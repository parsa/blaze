//=================================================================================================
/*!
//  \file blaze/math/blas/Level2.h
//  \brief Header file for BLAS level 2 functions
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

#ifndef _BLAZE_MATH_BLAS_LEVEL2_H_
#define _BLAZE_MATH_BLAS_LEVEL2_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/BlasCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/ConstDataAccess.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/system/BLAS.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>


namespace blaze {

//=================================================================================================
//
//  BLAS LEVEL 2 FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS level 2 functions */
//@{
#if BLAZE_BLAS_MODE

BLAZE_ALWAYS_INLINE void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, int M, int N,
                               float alpha, const float* A, int lda, const float* x, int incX,
                               float beta, float* y, int incY );

BLAZE_ALWAYS_INLINE void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, int M, int N,
                               double alpha, const double* A, int lda, const double* x, int incX,
                               double beta, double* y, int incY );

BLAZE_ALWAYS_INLINE void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, int M, int N,
                               complex<float> alpha, const complex<float>* A, int lda,
                               const complex<float>* x, int incX, complex<float> beta,
                               complex<float>* y, int incY );

BLAZE_ALWAYS_INLINE void gemv( CBLAS_ORDER layout, CBLAS_TRANSPOSE transA, int M, int N,
                               complex<double> alpha, const complex<double>* A, int lda,
                               const complex<double>* x, int incX, complex<double> beta,
                               complex<double>* y, int incY );

template< typename VT1, typename MT1, bool SO, typename VT2, typename ST >
BLAZE_ALWAYS_INLINE void gemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                               const DenseVector<VT2,false>& x, ST alpha, ST beta );

template< typename VT1, typename VT2, typename MT1, bool SO, typename ST >
BLAZE_ALWAYS_INLINE void gemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                               const DenseMatrix<MT1,SO>& A, ST alpha, ST beta );

BLAZE_ALWAYS_INLINE void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                               CBLAS_DIAG diag, int N, const float* A, int lda, float* X,
                               int incX );

BLAZE_ALWAYS_INLINE void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                               CBLAS_DIAG diag, int N, const double* A, int lda, double* X,
                               int incX );

BLAZE_ALWAYS_INLINE void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                               CBLAS_DIAG diag, int N, const complex<float>* A, int lda,
                               complex<float>* X, int incX );

BLAZE_ALWAYS_INLINE void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                               CBLAS_DIAG diag, int N, const complex<double>* A, int lda,
                               complex<double>* X, int incX );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void trmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                               CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void trmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                               CBLAS_UPLO uplo );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for single precision operands
//        (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup math
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param M The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param N The number of columns of matrix \a A \f$[0..\infty)\f$.
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
*/
BLAZE_ALWAYS_INLINE void gemv( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, int M, int N,
                               float alpha, const float* A, int lda, const float* x, int incX,
                               float beta, float* y, int incY )
{
   cblas_sgemv( order, transA, M, N, alpha, A, lda, x, incX, beta, y, incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for double precision operands
//        (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup math
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param M The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param N The number of columns of matrix \a A \f$[0..\infty)\f$.
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
*/
BLAZE_ALWAYS_INLINE void gemv( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, int M, int N,
                               double alpha, const double* A, int lda, const double* x, int incX,
                               double beta, double* y, int incY )
{
   cblas_dgemv( order, transA, M, N, alpha, A, lda, x, incX, beta, y, incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for single precision complex
//        operands (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup math
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param M The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param N The number of columns of matrix \a A \f$[0..\infty)\f$.
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
*/
BLAZE_ALWAYS_INLINE void gemv( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, int M, int N,
                               complex<float> alpha, const complex<float>* A, int lda,
                               const complex<float>* x, int incX, complex<float> beta,
                               complex<float>* y, int incY )
{
   cblas_cgemv( order, transA, M, N, &alpha, A, lda, x, incX, &beta, y, incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for double precision complex
//        operands (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup math
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param M The number of rows of matrix \a A \f$[0..\infty)\f$.
// \param N The number of columns of matrix \a A \f$[0..\infty)\f$.
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
*/
BLAZE_ALWAYS_INLINE void gemv( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, int M, int N,
                               complex<double> alpha, const complex<double>* A, int lda,
                               const complex<double>* x, int incX, complex<double> beta,
                               complex<double>* y, int incY )
{
   cblas_zgemv( order, transA, M, N, &alpha, A, lda, x, incX, &beta, y, incY );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication
//        (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The left-hand side dense matrix operand.
// \param x The right-hand side dense vector operand.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \return void
//
// This function performs the dense matrix/dense vector multiplication based on the BLAS gemv()
// functions. Note that the function only works for vectors and matrices with \c float, \c double,
// \c complex<float>, or \c complex<double> element type. The attempt to call the function with
// vectors and matrices of any other element type results in a compile time error.
*/
template< typename VT1   // Type of the left-hand side target vector
        , typename MT1   // Type of the left-hand side matrix operand
        , bool SO        // Storage order of the left-hand side matrix operand
        , typename VT2   // Type of the right-hand side vector operand
        , typename ST >  // Type of the scalar factors
BLAZE_ALWAYS_INLINE void gemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                               const DenseVector<VT2,false>& x, ST alpha, ST beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename VT2::ElementType );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   gemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasNoTrans, M, N, alpha,
         (~A).data(), lda, (~x).data(), 1, beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/dense matrix multiplication
//        (\f$ \vec{y}^T=\alpha*\vec{x}^T*A+\beta*\vec{y}^T \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param x The left-hand side dense vector operand.
// \param A The right-hand side dense matrix operand.
// \param alpha The scaling factor for \f$ \vec{x}^T*A \f$.
// \param beta The scaling factor for \f$ \vec{y}^T \f$.
// \return void
//
// This function performs the transpose dense vector/dense matrix multiplication based on the
// BLAS gemv() functions. Note that the function only works for vectors and matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with vectors and matrices of any other element type results in a compile time error.
*/
template< typename VT1   // Type of the left-hand side target vector
        , typename VT2   // Type of the left-hand side vector operand
        , typename MT1   // Type of the right-hand side matrix operand
        , bool SO        // Storage order of the right-hand side matrix operand
        , typename ST >  // Type of the scalar factors
BLAZE_ALWAYS_INLINE void gemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                               const DenseMatrix<MT1,SO>& A, ST alpha, ST beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename VT2::ElementType );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   gemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasTrans, M, N, alpha,
         (~A).data(), lda, (~x).data(), 1, beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication for single
//        precision operands (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param diag Specifies whether \a A is unitriangular (\a CblasNonUnit or \a CblasUnit).
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param a Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \return void
//
// This function performs the multiplication of a single precision triangular matrix by a vector
// based on the cblas_strmv() function.
*/
BLAZE_ALWAYS_INLINE void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                               CBLAS_DIAG diag, int n, const float* a, int lda, float* x,
                               int incX )
{
   cblas_strmv( order, uplo, transA, diag, n, a, lda, x, incX );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication for double
//        precision operands (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param diag Specifies whether \a A is unitriangular (\a CblasNonUnit or \a CblasUnit).
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param a Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \return void
//
// This function performs the multiplication of a double precision triangular matrix by a vector
// based on the cblas_dtrmv() function.
*/
BLAZE_ALWAYS_INLINE void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                               CBLAS_DIAG diag, int n, const double* a, int lda, double* x,
                               int incX )
{
   cblas_dtrmv( order, uplo, transA, diag, n, a, lda, x, incX );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication for single
//        precision complex operands (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param diag Specifies whether \a A is unitriangular (\a CblasNonUnit or \a CblasUnit).
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param a Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \return void
//
// This function performs the multiplication of a single precision complex triangular matrix by a
// vector based on the cblas_ctrmv() function.
*/
BLAZE_ALWAYS_INLINE void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                               CBLAS_DIAG diag, int n, const complex<float>* a, int lda,
                               complex<float>* x, int incX )
{
   cblas_ctrmv( order, uplo, transA, diag, n, a, lda, x, incX );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication for double
//        precision complex operands (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param diag Specifies whether \a A is unitriangular (\a CblasNonUnit or \a CblasUnit).
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param a Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param x Pointer to the first element of vector \a x.
// \param incX The stride within vector \a x.
// \return void
//
// This function performs the multiplication of a double precision complex triangular matrix by a
// vector based on the cblas_ztrmv() function.
*/
BLAZE_ALWAYS_INLINE void trmv( CBLAS_ORDER order, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                               CBLAS_DIAG diag, int n, const complex<double>* a, int lda,
                               complex<double>* x, int incX )
{
   cblas_ztrmv( order, uplo, transA, diag, n, a, lda, x, incX );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication
//        (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a triangular matrix by a vector based on the BLAS
// trmv() functions. Note that the function only works for vectors and matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with vectors and matrices of any other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void trmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                               CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   trmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
         uplo, CblasNoTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/triangular dense matrix multiplication
//        (\f$ \vec{y}^T=\vec{y}^T*A \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a vector and a triangular matrix based on the BLAS
// trmv() functions. Note that the function only works for vectors and matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with vectors and matrices of any other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void trmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                               CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   trmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
         uplo, CblasTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
