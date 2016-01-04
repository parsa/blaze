//=================================================================================================
/*!
//  \file blaze/math/lapack/QR.h
//  \brief Header file for LAPACK QR decomposition functions
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

#ifndef _BLAZE_MATH_LAPACK_QR_H_
#define _BLAZE_MATH_LAPACK_QR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BlasCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Complex.h>
#include <blaze/util/constraints/Double.h>
#include <blaze/util/constraints/Float.h>
#include <blaze/util/StaticAssert.h>
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

void sgeqrf_( int* m, int* n, float*  a, int* lda, float*  tau, float*  work, int* lwork, int* info );
void dgeqrf_( int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info );
void cgeqrf_( int* m, int* n, float*  a, int* lda, float*  tau, float*  work, int* lwork, int* info );
void zgeqrf_( int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info );

}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK QR DECOMPOSITION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK QR decomposition functions */
//@{
void geqrf( int* m, int* n, float* a, int* lda, float* tau, float*  work, int* lwork, int* info );

void geqrf( int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info );

void geqrf( int* m, int* n, complex<float>* a, int* lda, complex<float>* tau,
            complex<float>* work, int* lwork, int* info );

void geqrf( int* m, int* n, complex<double>* a, int* lda, complex<double>* tau,
            complex<double>* work, int* lwork, int* info );

template< typename MT, bool SO >
inline void geqrf( DenseMatrix<MT,SO>& A, typename MT::ElementType* tau );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QR decomposition of the given dense single precision matrix.
// \ingroup lapack
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( M, N ).
// \param work Auxiliary array; size >= max( 1, lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix QR decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK sgeqrf() function. The decomposition has the form

                              \f[ A = Q \dot R, \f]\n

// where the \c Q is represented as a product of elementary reflectors

                  \f[ Q = H(1) H(2) . . . H(k), with k = min(m,n).\f]\n

// Each H(i) has the form

                          \f[ H(i) = I - tau * v * v^T, \f]\n

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(0:i-1) = 0</tt> and
// <tt>v(i) = 1</tt>. <tt>v(i+1:m)</tt> is stored on exit in <tt>A(i+1:m,i)</tt>, and \c tau
// in \c tau(i). Thus on exit the elements on and above the diagonal of the matrix contain the
// min(M,N)-by-N upper trapezoidal matrix R (R is upper triangular if m >= n); the elements
// below the diagonal, with the array \a tau, represent the orthogonal matrix Q as a product
// of min(M,N) elementary reflectors.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the sgeqrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \a A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void geqrf( int* m, int* n, float* a, int* lda, float* tau, float* work, int* lwork, int* info )
{
   sgeqrf_( m, n, a, lda, tau, work, lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QR decomposition of the given dense double precision matrix.
// \ingroup lapack
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( M, N ).
// \param work Auxiliary array; size >= max( 1, lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix QR decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK dgeqrf() function. The decomposition has the form

                              \f[ A = Q \dot R, \f]\n

// where the \c Q is represented as a product of elementary reflectors

                  \f[ Q = H(1) H(2) . . . H(k), with k = min(m,n).\f]\n

// Each H(i) has the form

                          \f[ H(i) = I - tau * v * v^T, \f]\n

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(0:i-1) = 0</tt> and
// <tt>v(i) = 1</tt>. <tt>v(i+1:m)</tt> is stored on exit in <tt>A(i+1:m,i)</tt>, and \c tau
// in \c tau(i). Thus on exit the elements on and above the diagonal of the matrix contain the
// min(M,N)-by-N upper trapezoidal matrix R (R is upper triangular if m >= n); the elements
// below the diagonal, with the array \a tau, represent the orthogonal matrix Q as a product
// of min(M,N) elementary reflectors.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the dgeqrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \a A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void geqrf( int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info )
{
   dgeqrf_( m, n, a, lda, tau, work, lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QR decomposition of the given dense single precision complex matrix.
// \ingroup lapack
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( M, N ).
// \param work Auxiliary array; size >= max( 1, lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix QR decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK cgeqrf() function. The decomposition has the form

                              \f[ A = Q \dot R, \f]\n

// where the \c Q is represented as a product of elementary reflectors

                  \f[ Q = H(1) H(2) . . . H(k), with k = min(m,n).\f]\n

// Each H(i) has the form

                          \f[ H(i) = I - tau * v * v^T, \f]\n

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(0:i-1) = 0</tt> and
// <tt>v(i) = 1</tt>. <tt>v(i+1:m)</tt> is stored on exit in <tt>A(i+1:m,i)</tt>, and \c tau
// in \c tau(i). Thus on exit the elements on and above the diagonal of the matrix contain the
// min(M,N)-by-N upper trapezoidal matrix R (R is upper triangular if m >= n); the elements
// below the diagonal, with the array \a tau, represent the orthogonal matrix Q as a product
// of min(M,N) elementary reflectors.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the cgeqrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \a A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void geqrf( int* m, int* n, complex<float>* a, int* lda, complex<float>* tau,
                   complex<float>* work, int* lwork, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cgeqrf_( m, n, reinterpret_cast<float*>( a ), lda, reinterpret_cast<float*>( tau ),
            reinterpret_cast<float*>( work ), lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QR decomposition of the given dense double precision complex matrix.
// \ingroup lapack
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( M, N ).
// \param work Auxiliary array; size >= max( 1, lwork ).
// \param lwork The dimension of the array \a work; size >= max( 1, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix QR decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK zgeqrf() function. The decomposition has the form

                              \f[ A = Q \dot R, \f]\n

// where the \c Q is represented as a product of elementary reflectors

                  \f[ Q = H(1) H(2) . . . H(k), with k = min(m,n).\f]\n

// Each H(i) has the form

                          \f[ H(i) = I - tau * v * v^T, \f]\n

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(0:i-1) = 0</tt> and
// <tt>v(i) = 1</tt>. <tt>v(i+1:m)</tt> is stored on exit in <tt>A(i+1:m,i)</tt>, and \c tau
// in \c tau(i). Thus on exit the elements on and above the diagonal of the matrix contain the
// min(M,N)-by-N upper trapezoidal matrix R (R is upper triangular if m >= n); the elements
// below the diagonal, with the array \a tau, represent the orthogonal matrix Q as a product
// of min(M,N) elementary reflectors.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: The i-th argument had an illegal value.
//
// For more information on the zgeqrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \a A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void geqrf( int* m, int* n, complex<double>* a, int* lda, complex<double>* tau,
                   complex<double>* work, int* lwork, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   zgeqrf_( m, n, reinterpret_cast<double*>( a ), lda, reinterpret_cast<double*>( tau ),
            reinterpret_cast<double*>( work ), lwork, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the QR decomposition of the given dense matrix.
// \ingroup lapack
//
// \param A The matrix to be decomposed.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( M, N ).
// \return void
//
// This function performs the dense matrix QR decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK geqrf() functions. Note that this function can only be used for general,
// non-adapted matrices with \c float, \c double, \c complex<float>, or \c complex<double> element
// type. The attempt to call the function with any adapted matrix or matrices of any other element
// type results in a compile time error!\n
//
// The decomposition has the form

                              \f[ A = Q \dot R, \f]\n

// where the \c Q is represented as a product of elementary reflectors

                  \f[ Q = H(1) H(2) . . . H(k), with k = min(m,n).\f]\n

// Each H(i) has the form

                          \f[ H(i) = I - tau * v * v^T, \f]\n

// where \c tau is a real scalar, and \c v is a real vector with <tt>v(0:i-1) = 0</tt> and
// <tt>v(i) = 1</tt>. <tt>v(i+1:m)</tt> is stored on exit in <tt>A(i+1:m,i)</tt>, and \c tau
// in \c tau(i). Thus on exit the elements on and above the diagonal of the matrix contain the
// min(M,N)-by-N upper trapezoidal matrix R (R is upper triangular if m >= n); the elements
// below the diagonal, with the array \a tau, represent the orthogonal matrix Q as a product
// of min(M,N) elementary reflectors.
//
// For more information on the geqrf() functions (i.e. sgeqrf(), dgeqrf(), cgeqrf(), and zgeqrf())
// see the LAPACK online documentation browser:
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
inline void geqrf( DenseMatrix<MT,SO>& A, typename MT::ElementType* tau )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

   typedef typename MT::ElementType  ET;

   int m   ( boost::numeric_cast<int>( (~A).rows()    ) );
   int n   ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda ( boost::numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   if( m == 0 || n == 0 ) {
      return;
   }

   int lwork( n*lda );
   const UniqueArray<ET> work( new ET[lwork] );

   geqrf( &m, &n, (~A).data(), &lda, tau, work.get(), &lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for QR decomposition" );
}
//*************************************************************************************************

} // namespace blaze

#endif
