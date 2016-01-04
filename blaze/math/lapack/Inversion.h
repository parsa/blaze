//=================================================================================================
/*!
//  \file blaze/math/lapack/Inversion.h
//  \brief Header file for LAPACK matrix inversion functions
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

#ifndef _BLAZE_MATH_LAPACK_INVERSION_H_
#define _BLAZE_MATH_LAPACK_INVERSION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BlasCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/lapack/PLU.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/Complex.h>
#include <blaze/util/constraints/Double.h>
#include <blaze/util/constraints/Float.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/typetraits/IsComplexDouble.h>
#include <blaze/util/typetraits/IsComplexFloat.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsFloat.h>
#include <blaze/util/UniqueArray.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
extern "C" {

void sgetri_( int* n, float*  a, int* lda, int* ipiv, float*  work, int* lwork, int* info );
void dgetri_( int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info );
void cgetri_( int* n, float*  a, int* lda, int* ipiv, float*  work, int* lwork, int* info );
void zgetri_( int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info );

void spotri_( char* uplo, int* n, float*  a, int* lda, int* info );
void dpotri_( char* uplo, int* n, double* a, int* lda, int* info );
void cpotri_( char* uplo, int* n, float*  a, int* lda, int* info );
void zpotri_( char* uplo, int* n, double* a, int* lda, int* info );

}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK LU-BASED INVERSION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK LU-based inversion functions */
//@{
inline void getri( int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info );

inline void getri( int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info );

inline void getri( int* n, complex<float>* a, int* lda, int* ipiv,
                   complex<float>* work, int* lwork, int* info );

inline void getri( int* n, complex<double>* a, int* lda, int* ipiv,
                   complex<double>* work, int* lwork, int* info );

template< typename MT, bool SO >
inline void getri( DenseMatrix<MT,SO>& A, const int* ipiv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense single precision matrix.
// \ingroup lapack
//
// \param n The number of rows/columns of the matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix to be inverted.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \param work Auxiliary array; size >= max( 1, lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK sgetri() function for
// single precision matrices that have already been factorized by the sgetrf() function. The
// \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the inversion could not be computed since U(i,i) is exactly zero.
//
// If the function exits successfully (i.e. \a info = 0) then the first element of the \a work
// array returns the optimal \a lwork. For optimal performance \a lwork >= N*NB, where NB is the
// optimal blocksize returned by the LAPACK function ilaenv(). If \a lwork = -1 then a workspace
// query is assumed. The function only calculates the optimal size of the \a work array and
// returns this value as the first entry of the \a work array.
//
// For more information on the sgetri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void getri( int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info )
{
   sgetri_( n, a, lda, ipiv, work, lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense double precision matrix.
// \ingroup lapack
//
// \param n The number of rows/columns of the matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix to be inverted.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \param work Auxiliary array; size >= max( 1, lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK dgetri() function for
// double precision matrices that have already been factorized by the dgetrf() function. The
// \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the inversion could not be computed since U(i,i) is exactly zero.
//
// If the function exits successfully (i.e. \a info = 0) then the first element of the \a work
// array returns the optimal \a lwork. For optimal performance \a lwork >= N*NB, where NB is the
// optimal blocksize returned by the LAPACK function ilaenv(). If \a lwork = -1 then a workspace
// query is assumed. The function only calculates the optimal size of the \a work array and
// returns this value as the first entry of the \a work array.
//
// For more information on the sgetri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void getri( int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info )
{
   dgetri_( n, a, lda, ipiv, work, lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense single precision complex matrix.
// \ingroup lapack
//
// \param n The number of rows/columns of the matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix to be inverted.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \param work Auxiliary array; size >= max( 1, lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK cgetri() function for
// single precision complex matrices that have already been factorized by the cgetrf() function.
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the inversion could not be computed since U(i,i) is exactly zero.
//
// If the function exits successfully (i.e. \a info = 0) then the first element of the \a work
// array returns the optimal \a lwork. For optimal performance \a lwork >= N*NB, where NB is the
// optimal blocksize returned by the LAPACK function ilaenv(). If \a lwork = -1 then a workspace
// query is assumed. The function only calculates the optimal size of the \a work array and
// returns this value as the first entry of the \a work array.
//
// For more information on the sgetri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void getri( int* n, complex<float>* a, int* lda, int* ipiv,
                   complex<float>* work, int* lwork, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cgetri_( n, reinterpret_cast<float*>( a ), lda, ipiv,
            reinterpret_cast<float*>( work ), lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense double precision complex matrix.
// \ingroup lapack
//
// \param n The number of rows/columns of the matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix to be inverted.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \param work Auxiliary array; size >= max( 1, lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK zgetri() function for
// double precision complex matrices that have already been factorized by the zgetrf() function.
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the inversion could not be computed since U(i,i) is exactly zero.
//
// If the function exits successfully (i.e. \a info = 0) then the first element of the \a work
// array returns the optimal \a lwork. For optimal performance \a lwork >= N*NB, where NB is the
// optimal blocksize returned by the LAPACK function ilaenv(). If \a lwork = -1 then a workspace
// query is assumed. The function only calculates the optimal size of the \a work array and
// returns this value as the first entry of the \a work array.
//
// For more information on the sgetri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void getri( int* n, complex<double>* a, int* lda, int* ipiv,
                   complex<double>* work, int* lwork, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   zgetri_( n, reinterpret_cast<double*>( a ), lda, ipiv,
            reinterpret_cast<double*>( work ), lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense matrix.
// \ingroup lapack
//
// \param A The matrix to be inverted.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function performs the dense matrix inversion based on the LAPACK getri() functions for
// matrices that have already been factorized by the getrf() functions. Note that the function
// only works for general, non-adapted matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with adaptors or matrices
// of any other element type results in a compile time error!
//
// The function fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the getri() functions (i.e. sgetri(), dgetri(), cgetri(), and zgetri())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \c A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void getri( DenseMatrix<MT,SO>& A, const int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

   typedef typename MT::ElementType  ET;

   if( !isSquare( ~A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   int n   ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda ( boost::numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   if( n == 0 )
      return;

   int lwork( n*lda );
   const UniqueArray<ET> work( new ET[lwork] );

   getri( &n, (~A).data(), &lda, const_cast<int*>( ipiv ), work.get(), &lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for matrix inversion" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK CHOLESKY-BASED INVERSION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK Cholesky-based inversion functions */
//@{
inline void potri( char* uplo, int* n, float*  a, int* lda, int* info );

inline void potri( char* uplo, int* n, double* a, int* lda, int* info );

inline void potri( char* uplo, int* n, complex<float>*  a, int* lda, int* info );

inline void potri( char* uplo, int* n, complex<double>* a, int* lda, int* info );

template< typename MT, bool SO >
inline void potri( DenseMatrix<MT,SO>& A, char uplo );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense positive-definite single precision
//        matrix.
// \ingroup lapack
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of the matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK spotri() function for
// positive-definite single precision matrices that have already been factorized by the spotrf()
// function. The resulting symmetric inverse of \a A is stored either in the lower part of \a A
// (\a uplo = \c 'L') or in the upper part (\a uplo = \c 'U').
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If \a info = -i, the i-th argument had an illegal value.
//   - > 0: If \a info = i, element (i,i) of U or L is zero and the inverse could not be computed.
//
// For more information on the spotri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void potri( char* uplo, int* n, float* a, int* lda, int* info )
{
   spotri_( uplo, n, a, lda, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense positive-definite double precision
//        matrix.
// \ingroup lapack
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of the matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK dpotri() function for
// positive-definite double precision matrices that have already been factorized by the dpotrf()
// function. The resulting symmetric inverse of \a A is stored either in the lower part of \a A
// (\a uplo = \c 'L') or in the upper part (\a uplo = \c 'U').
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If \a info = -i, the i-th argument had an illegal value.
//   - > 0: If \a info = i, element (i,i) of U or L is zero and the inverse could not be computed.
//
// For more information on the spotri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void potri( char* uplo, int* n, double* a, int* lda, int* info )
{
   dpotri_( uplo, n, a, lda, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense positive-definite single precision
//        complex matrix.
// \ingroup lapack
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of the matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK cpotri() function for
// positive-definite single precision complex matrices that have already been factorized by the
// cpotrf() function. The resulting symmetric inverse of \a A is stored either in the lower part
// of \a A (\a uplo = \c 'L') or in the upper part (\a uplo = \c 'U').
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If \a info = -i, the i-th argument had an illegal value.
//   - > 0: If \a info = i, element (i,i) of U or L is zero and the inverse could not be computed.
//
// For more information on the cpotri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void potri( char* uplo, int* n, complex<float>* a, int* lda, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cpotri_( uplo, n, reinterpret_cast<float*>( a ), lda, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense positive-definite double precision
//        complex matrix.
// \ingroup lapack
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of the matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK zpotri() function for
// positive-definite double precision complex matrices that have already been factorized by the
// zpotrf() function. The resulting symmetric inverse of \a A is stored either in the lower part
// of \a A (\a uplo = \c 'L') or in the upper part (\a uplo = \c 'U').
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If \a info = -i, the i-th argument had an illegal value.
//   - > 0: If \a info = i, element (i,i) of U or L is zero and the inverse could not be computed.
//
// For more information on the zpotri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void potri( char* uplo, int* n, complex<double>* a, int* lda, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   zpotri_( uplo, n, reinterpret_cast<double*>( a ), lda, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense positive-definite matrix.
// \ingroup lapack
//
// \param A The positive-definite matrix to be inverted.
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \return void
// \exception std::invalid_argument Invalid argument provided.
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function performs the dense matrix inversion based on the LAPACK potri() functions for
// positive-definite matrices that have already been factorized by the potrf() functions. The
// resulting symmetric inverse of the given matrix \a A is stored either in the lower part of
// \a A (\a uplo = \c 'L') or in the upper part (\a uplo = \c 'U'). Note that the function only
// works for general, non-adapted matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with adaptors or matrices
// of any other element type results in a compile time error!
//
// The function fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given \a uplo argument is neither 'L' nor 'U';
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the potri() functions (i.e. spotri(), dpotri(), cpotri(), and zpotri())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \c A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void potri( DenseMatrix<MT,SO>& A, char uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

   if( !isSquare( ~A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   if( uplo != 'L' && uplo != 'U' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid uplo argument provided" );
   }

   int n   ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda ( boost::numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   if( n == 0 ) {
      return;
   }

   if( IsRowMajorMatrix<MT>::value ) {
      ( uplo == 'L' )?( uplo = 'U' ):( uplo = 'L' );
   }

   potri( &uplo, &n, (~A).data(), &lda, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for matrix inversion" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
