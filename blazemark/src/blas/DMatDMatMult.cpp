//=================================================================================================
/*!
//  \file src/blas/DMatDMatMult.cpp
//  \brief Source file for the BLAS dense matrix/dense matrix multiplication kernel
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
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/blas/DMatDMatMult.h>
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
/*!\brief Kernel function for single-precision matrices.
//
// \param Order Whether matrices are row major order (C-Style) for column major order (Fortran-style).
// \param TransA Whether to transpose matrix A.
// \param TransB Whether to transpose matrix B.
// \param M Rows in matrices A and C.
// \param N Columns in Matrices B and C.
// \param K Columns in matrix A and Rows in matrix B.
// \param alpha Scalar factor for \f$ op(A)op(B) \f$.
// \param A Pointer to the first element of matrix A.
// \param lda The size of the first dimension of matrix A.
// \param B Pointer to the first element of matrix B.
// \param ldb The size of the first dimension of matrix B.
// \param beta Scalar factor for \f$ C \f$.
// \param C Pointer to the first element of matrix C.
// \param ldc The size of the first dimension of matrix C.
// \return void
*/
inline void gemm( const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                  const CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
                  const float alpha, const float *A, const int lda, const float *B, const int ldb,
                  const float beta, float *C, const int ldc )
{
   cblas_sgemm( Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Kernel function for double-precision matrices.
//
// \param Order Whether matrices are row major order (C-Style) for column major order (Fortran-style).
// \param TransA Whether to transpose matrix A.
// \param TransB Whether to transpose matrix B.
// \param M Rows in matrices A and C.
// \param N Columns in Matrices B and C.
// \param K Columns in matrix A and Rows in matrix B.
// \param alpha Scalar factor for \f$ op(A)op(B) \f$.
// \param A Pointer to the first element of matrix A.
// \param lda The size of the first dimension of matrix A.
// \param B Pointer to the first element of matrix B.
// \param ldb The size of the first dimension of matrix B.
// \param beta Scalar factor for \f$ C \f$.
// \param C Pointer to the first element of matrix C.
// \param ldc The size of the first dimension of matrix C.
// \return void
*/
inline void gemm( const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                  const CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
                  const double alpha, const double *A, const int lda, const double *B, const int ldb,
                  const double beta, double *C, const int ldc )
{
   cblas_dgemm( Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc );
}
//*************************************************************************************************




//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief BLAS dense matrix/dense matrix multiplication kernel.
//
// \param N The number of rows and columns of the matrices.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the dense matrix/dense matrix multiplication by means of
// BLAS functionality.
*/
double dmatdmatmult( size_t N, size_t steps )
{
   using ::blazemark::real;
   using ::blaze::rowMajor;

   ::blaze::setSeed( seed );

   ::blaze::DynamicMatrix<real,rowMajor> A( N, N ), B( N, N ), C( N, N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         A(i,j) = ::blaze::rand<real>();
         B(i,j) = ::blaze::rand<real>();
      }
   }

   gemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, real(1),
         A.data(), A.spacing(), B.data(), B.spacing(), real(0), C.data(), C.spacing() );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         gemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, real(1),
               A.data(), A.spacing(), B.data(), B.spacing(), real(0), C.data(), C.spacing() );
      }
      timer.end();

      if( C.rows() != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " BLAS kernel 'dmatdmatmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blas

} // namespace blazemark
