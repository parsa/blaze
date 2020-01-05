//=================================================================================================
/*!
//  \file blaze/math/blas/cblas/gemv.h
//  \brief Header file for the CBLAS gemv wrapper functions
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_BLAS_CBLAS_GEMV_H_
#define _BLAZE_MATH_BLAS_CBLAS_GEMV_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/Types.h>
#include <blaze/system/BLAS.h>
#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>


//=================================================================================================
//
//  BLAS FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if !defined(INTEL_MKL_VERSION)
extern "C" {

void sgemv_( char* trans, blaze::blas_int_t* m, blaze::blas_int_t* n, float* alpha,
             float* A, blaze::blas_int_t* lda, float* x, blaze::blas_int_t* incX,
             float* beta, float* y, blaze::blas_int_t* incY,
             blaze::fortran_charlen_t ntransA, blaze::fortran_charlen_t ntransB );
void dgemv_( char* trans, blaze::blas_int_t* m, blaze::blas_int_t* n, double* alpha,
             double* A, blaze::blas_int_t* lda, double* x, blaze::blas_int_t* incX,
             double* beta, double* y, blaze::blas_int_t* incY,
             blaze::fortran_charlen_t ntransA, blaze::fortran_charlen_t ntransB );
void cgemv_( char* trans, blaze::blas_int_t* m, blaze::blas_int_t* n, float* alpha,
             float* A, blaze::blas_int_t* lda, float* x, blaze::blas_int_t* incX,
             float* beta, float* y, blaze::blas_int_t* incY,
             blaze::fortran_charlen_t ntransA, blaze::fortran_charlen_t ntransB );
void zgemv_( char* trans, blaze::blas_int_t* m, blaze::blas_int_t* n, double* alpha,
             double* A, blaze::blas_int_t* lda, double* x, blaze::blas_int_t* incX,
             double* beta, double* y, blaze::blas_int_t* incY,
             blaze::fortran_charlen_t ntransA, blaze::fortran_charlen_t ntransB );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  BLAS GENERAL MATRIX/VECTOR MULTIPLICATION FUNCTIONS (GEMV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS general matrix/vector multiplication functions (gemv) */
//@{
void gemv( char trans, blas_int_t m, blas_int_t n, float alpha, const float* A,
           blas_int_t lda, const float* x, blas_int_t incX, float beta, float* y,
           blas_int_t incY );

void gemv( char trans, blas_int_t m, blas_int_t n, double alpha, const double* A,
           blas_int_t lda, const double* x, blas_int_t incX, double beta, double* y,
           blas_int_t incY );

void gemv( char trans, blas_int_t m, blas_int_t n, complex<float> alpha,
           const complex<float>* A, blas_int_t lda, const complex<float>* x,
           blas_int_t incX, complex<float> beta, complex<float>* y, blas_int_t incY );

void gemv( char trans, blas_int_t m, blas_int_t n, complex<double> alpha,
           const complex<double>* A, blas_int_t lda, const complex<double>* x,
           blas_int_t incX, complex<double> beta, complex<double>* y, blas_int_t incY );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication with single precision
//        column-major matrices (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param trans \c 'N' to use \a A, \c 'T' or \c 'C' to use \a A^T.
// \param m The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a A \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for general single precision
// matrices based on the BLAS sgemv() function (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
//
// For more information on the sgemv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemv( char trans, blas_int_t m, blas_int_t n, float alpha, const float* A,
                  blas_int_t lda, const float* x, blas_int_t incX, float beta, float* y,
                  blas_int_t incY )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   sgemv_( &trans, &m, &n, &alpha, const_cast<float*>( A ), &lda,
           const_cast<float*>( x ), &incX, &beta, y, &incY
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication with double precision
//        column-major matrices (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param trans \c 'N' to use \a A, \c 'T' or \c 'C' to use \a A^T.
// \param m The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a A \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for general double precision
// matrices based on the BLAS dgemv() function (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
//
// For more information on the dgemv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemv( char trans, blas_int_t m, blas_int_t n, double alpha, const double* A,
                  blas_int_t lda, const double* x, blas_int_t incX, double beta, double* y,
                  blas_int_t incY )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   dgemv_( &trans, &m, &n, &alpha, const_cast<double*>( A ), &lda,
           const_cast<double*>( x ), &incX, &beta, y, &incY
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication with single precision
//        complex column-major matrices (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param trans \c 'N' to use \a A, \c 'T' or \c 'C' to use \a A^T.
// \param m The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a A \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for general single
// precision complex matrices based on the BLAS cgemv() function
// (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
//
// For more information on the cgemv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemv( char trans, blas_int_t m, blas_int_t n, complex<float> alpha,
                  const complex<float>* A, blas_int_t lda, const complex<float>* x,
                  blas_int_t incX, complex<float> beta, complex<float>* y, blas_int_t incY )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   cgemv_( &trans, &m, &n, reinterpret_cast<ET*>( &alpha ),
           const_cast<ET*>( reinterpret_cast<const ET*>( A ) ), &lda,
           const_cast<ET*>( reinterpret_cast<const ET*>( x ) ), &incX,
           reinterpret_cast<ET*>( &beta ), reinterpret_cast<ET*>( y ), &incY
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication with double precision
//        complex column-major matrices (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param trans \c 'N' to use \a A, \c 'T' or \c 'C' to use \a A^T.
// \param m The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a A \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for general double
// precision complex matrices based on the BLAS zgemv() function
// (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
//
// For more information on the zgemv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemv( char trans, blas_int_t m, blas_int_t n, complex<double> alpha,
                  const complex<double>* A, blas_int_t lda, const complex<double>* x,
                  blas_int_t incX, complex<double> beta, complex<double>* y, blas_int_t incY )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zgemv_( &trans, &m, &n, reinterpret_cast<ET*>( &alpha ),
           const_cast<ET*>( reinterpret_cast<const ET*>( A ) ), &lda,
           const_cast<ET*>( reinterpret_cast<const ET*>( x ) ), &incX,
           reinterpret_cast<ET*>( &beta ), reinterpret_cast<ET*>( y ), &incY
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************




//=================================================================================================
//
//  CBLAS WRAPPER FUNCTIONS (GEMV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS wrapper functions (gemv) */
//@{
#if BLAZE_BLAS_MODE

void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
           float alpha, const float* A, blas_int_t lda, const float* x,
           blas_int_t incX, float beta, float* y, blas_int_t incY );

void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
           double alpha, const double* A, blas_int_t lda, const double* x,
           blas_int_t incX, double beta, double* y, blas_int_t incY );

void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
           complex<float> alpha, const complex<float>* A, blas_int_t lda,
           const complex<float>* x, blas_int_t incX, complex<float> beta,
           complex<float>* y, blas_int_t incY );

void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
           complex<double> alpha, const complex<double>* A, blas_int_t lda,
           const complex<double>* x, blas_int_t incX, complex<double> beta,
           complex<double>* y, blas_int_t incY );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for single precision operands
//        (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param m The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a A \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for single precision
// operands based on the BLAS cblas_sgemv() function.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemv( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
                  float alpha, const float* A, blas_int_t lda, const float* x,
                  blas_int_t incX, float beta, float* y, blas_int_t incY )
{
   cblas_sgemv( order, transA, m, n, alpha, A, lda, x, incX, beta, y, incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for double precision operands
//        (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param m The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a A \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for double precision
// operands based on the BLAS cblas_dgemv() function.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemv( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
                  double alpha, const double* A, blas_int_t lda, const double* x,
                  blas_int_t incX, double beta, double* y, blas_int_t incY )
{
   cblas_dgemv( order, transA, m, n, alpha, A, lda, x, incX, beta, y, incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for single precision complex
//        operands (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param m The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a A \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for single precision
// complex operands based on the BLAS cblas_cgemv() function.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemv( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
                  complex<float> alpha, const complex<float>* A, blas_int_t lda,
                  const complex<float>* x, blas_int_t incX, complex<float> beta,
                  complex<float>* y, blas_int_t incY )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cblas_cgemv( order, transA, m, n, reinterpret_cast<const float*>( &alpha ),
                reinterpret_cast<const float*>( A ), lda, reinterpret_cast<const float*>( x ),
                incX, reinterpret_cast<const float*>( &beta ), reinterpret_cast<float*>( y ), incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for double precision complex
//        operands (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param m The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a A \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for double precision
// complex operands based on the BLAS zblas_zgemv() function.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemv( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, blas_int_t m, blas_int_t n,
                  complex<double> alpha, const complex<double>* A, blas_int_t lda,
                  const complex<double>* x, blas_int_t incX, complex<double> beta,
                  complex<double>* y, blas_int_t incY )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   cblas_zgemv( order, transA, m, n, reinterpret_cast<const double*>( &alpha ),
                reinterpret_cast<const double*>( A ), lda, reinterpret_cast<const double*>( x ),
                incX, reinterpret_cast<const double*>( &beta ), reinterpret_cast<double*>( y ), incY );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
