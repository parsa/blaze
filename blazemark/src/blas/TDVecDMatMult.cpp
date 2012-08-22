//=================================================================================================
/*!
//  \file src/blas/TDVecDMatMult.cpp
//  \brief Source file for the BLAS transpose dense vector/dense matrix multiplication kernel
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
*/
//=================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iostream>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/blas/TDVecDMatMult.h>
#include <blazemark/system/BLAS.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Precision.h>


namespace blazemark {

namespace blas {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Kernel function for single-precision vectors and matrices.
//
// \param Order Whether the matrix is row major order (C-Style) for column major order (Fortran-style).
// \param TransA Whether to transpose matrix A.
// \param M Rows in matrices A.
// \param N Columns in matrix A.
// \param alpha Scalar factor for \f$ \alpha op(A)x \f$.
// \param A Pointer to the first element of matrix A.
// \param lda The size of the first dimension of matrix A.
// \param X Pointer to the first element of vector X.
// \param incX Use every incX'th element of vector X.
// \param beta Scalar factor for \f$ y \f$.
// \param Y Pointer to the first element of vector Y.
// \param incY Use every incY'th element of vector Y.
// \return void
*/
inline void gemv( const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                  const int M, const int N, const float alpha,
                  const float *A, const int lda, const float *X, const int incX,
                  const float beta, float *Y, const int incY )
{
   cblas_sgemv( Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Kernel function for double-precision vectors and matrices.
//
// \param Order Whether the matrix is row major order (C-Style) for column major order (Fortran-style).
// \param TransA Whether to transpose matrix A.
// \param M Rows in matrices A.
// \param N Columns in matrix A.
// \param alpha Scalar factor for \f$ \alpha op(A)x \f$.
// \param A Pointer to the first element of matrix A.
// \param lda The size of the first dimension of matrix A.
// \param X Pointer to the first element of vector X.
// \param incX Use every incX'th element of vector X.
// \param beta Scalar factor for \f$ y \f$.
// \param Y Pointer to the first element of vector Y.
// \param incY Use every incY'th element of vector Y.
// \return void
*/
inline void gemv( const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                  const int M, const int N, const double alpha,
                  const double *A, const int lda, const double *X, const int incX,
                  const double beta, double *Y, const int incY )
{
   cblas_dgemv( Order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY );
}
//*************************************************************************************************




//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief BLAS transpose dense vector/dense matrix multiplication kernel.
//
// \param N The number of rows and columns of the matrix and the size of the vector.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the transpose dense vector/dense matrix multiplication by
// means of BLAS functionality.
*/
double tdvecdmatmult( size_t N, size_t steps )
{
   using ::blazemark::real;
   using ::blaze::rowVector;
   using ::blaze::rowMajor;

   ::blaze::setSeed( seed );

   ::blaze::DynamicMatrix<real,rowMajor> A( N, N );
   ::blaze::DynamicVector<real,rowVector> a( N ), b( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         A(i,j) = ::blaze::rand<real>();
      }
   }

   for( size_t i=0UL; i<N; ++i ) {
      a[i] = ::blaze::rand<real>();
   }

   gemv( CblasRowMajor, CblasTrans, N, N, real(1),
         A.data(), A.spacing(), a.data(), 1, real(0), b.data(), 1 );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         gemv( CblasRowMajor, CblasTrans, N, N, real(1),
               A.data(), A.spacing(), a.data(), 1, real(0), b.data(), 1 );
      }
      timer.end();

      if( b.size() != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " BLAS kernel 'tdvecdmatmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blas

} // namespace blazemark
