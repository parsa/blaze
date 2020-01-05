//=================================================================================================
/*!
//  \file blaze/math/blas/cblas/axpy.h
//  \brief Header file for the CBLAS axpy wrapper functions
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

#ifndef _BLAZE_MATH_BLAS_CBLAS_AXPY_H_
#define _BLAZE_MATH_BLAS_CBLAS_AXPY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/Types.h>
#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>


//=================================================================================================
//
//  BLAS FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if !defined(INTEL_MKL_VERSION)
extern "C" {

void saxpy_( blaze::blas_int_t* n, float* alpha, float* x, blaze::blas_int_t* incX,
             float* y, blaze::blas_int_t* incY );

void daxpy_( blaze::blas_int_t* n, double* alpha, double* x, blaze::blas_int_t* incX,
             double* y, blaze::blas_int_t* incY );

void caxpy_( blaze::blas_int_t* n, float* alpha, float* x, blaze::blas_int_t* incX,
             float* y, blaze::blas_int_t* incY );

void zaxpy_( blaze::blas_int_t* n, double* alpha, double* x, blaze::blas_int_t* incX,
             double* y, blaze::blas_int_t* incY );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  BLAS SCALED VECTOR ADDITION (AXPY)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS scaled vector addition functions (axpy) */
//@{
void axpy( blas_int_t n, float alpha, const float* x, blas_int_t incX,
           float* y, blas_int_t incY );

void axpy( blas_int_t n, double alpha, const double* x, blas_int_t incX,
           double* y, blas_int_t incY );

void axpy( blas_int_t n, complex<float> alpha, const complex<float>* x, blas_int_t incX,
           complex<float>* y, blas_int_t incY );

void axpy( blas_int_t n, complex<double> alpha, const complex<double>* x, blas_int_t incX,
           complex<double>* y, blas_int_t incY );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief BLAS kernel for scaled dense vector addition for single precision operands
//        (\f$ \vec{y}+=\alpha*\vec{x} \f$).
// \ingroup blas
//
// \param n The size of the two dense vectors \a x and \a y \f$[0..\infty)\f$.
// \param alpha The scaling factor for the dense vector \a x.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the a scaled dense vector addition for single precision operands based
// on the BLAS saxpy() function (\f$ \vec{y}+=\alpha*\vec{x} \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void axpy( blas_int_t n, float alpha, const float* x, blas_int_t incX,
                  float* y, blas_int_t incY )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   saxpy_( &n, &alpha, const_cast<float*>( x ), &incX, y, &incY );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief BLAS kernel for scaled dense vector addition for double precision operands
//        (\f$ \vec{y}+=\alpha*\vec{x} \f$).
// \ingroup blas
//
// \param n The size of the two dense vectors \a x and \a y \f$[0..\infty)\f$.
// \param alpha The scaling factor for the dense vector \a x.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the a scaled dense vector addition for double precision operands based
// on the BLAS daxpy() function (\f$ \vec{y}+=\alpha*\vec{x} \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void axpy( blas_int_t n, double alpha, const double* x, blas_int_t incX,
                  double* y, blas_int_t incY )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   daxpy_( &n, &alpha, const_cast<double*>( x ), &incX, y, &incY );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief BLAS kernel for scaled dense vector addition for single precision complex operands
//        (\f$ \vec{y}+=\alpha*\vec{x} \f$).
// \ingroup blas
//
// \param n The size of the two dense vectors \a x and \a y \f$[0..\infty)\f$.
// \param alpha The scaling factor for the dense vector \a x.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the a scaled dense vector addition for single precision complex
// operands based on the BLAS caxpy() function (\f$ \vec{y}+=\alpha*\vec{x} \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void axpy( blas_int_t n, complex<float> alpha, const complex<float>* x, blas_int_t incX,
                  complex<float>* y, blas_int_t incY )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   caxpy_( &n, reinterpret_cast<float*>( &alpha ),
           const_cast<ET*>( reinterpret_cast<const ET*>( x ) ), &incX,
           reinterpret_cast<ET*>( y ), &incY );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief BLAS kernel for scaled dense vector addition for double precision complex operands
//        (\f$ \vec{y}+=\alpha*\vec{x} \f$).
// \ingroup blas
//
// \param n The size of the two dense vectors \a x and \a y \f$[0..\infty)\f$.
// \param alpha The scaling factor for the dense vector \a x.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \param y Pointer to the first element of vector \a y.
// \param incY The stride within vector \a y.
// \return void
//
// This function performs the a scaled dense vector addition for double precision complex
// operands based on the BLAS zaxpy() function (\f$ \vec{y}+=\alpha*\vec{x} \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void axpy( blas_int_t n, complex<double> alpha, const complex<double>* x, blas_int_t incX,
                  complex<double>* y, blas_int_t incY )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zaxpy_( &n, reinterpret_cast<double*>( &alpha ),
           const_cast<ET*>( reinterpret_cast<const ET*>( x ) ), &incX,
           reinterpret_cast<ET*>( y ), &incY );
}
//*************************************************************************************************

} // namespace blaze

#endif
