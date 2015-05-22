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
#include <blaze/util/constraints/Complex.h>
#include <blaze/util/constraints/Double.h>
#include <blaze/util/constraints/Float.h>


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

template< typename VT1, typename MT1, bool SO, typename VT2 >
BLAZE_ALWAYS_INLINE void sgemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                                const DenseVector<VT2,false>& x, float alpha, float beta );

template< typename VT1, typename VT2, typename MT1, bool SO >
BLAZE_ALWAYS_INLINE void sgemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                                const DenseMatrix<MT1,SO>& A, float alpha, float beta );

template< typename VT1, typename MT1, bool SO, typename VT2 >
BLAZE_ALWAYS_INLINE void dgemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                                const DenseVector<VT2,false>& x, double alpha, double beta );

template< typename VT1, typename VT2, typename MT1, bool SO >
BLAZE_ALWAYS_INLINE void dgemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                                const DenseMatrix<MT1,SO>& A, double alpha, double beta );

template< typename VT1, typename MT1, bool SO, typename VT2 >
BLAZE_ALWAYS_INLINE void cgemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                                const DenseVector<VT2,false>& x, complex<float> alpha, complex<float> beta );

template< typename VT1, typename VT2, typename MT1, bool SO >
BLAZE_ALWAYS_INLINE void cgemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                                const DenseMatrix<MT1,SO>& A, complex<float> alpha, complex<float> beta );

template< typename VT1, typename MT1, bool SO, typename VT2 >
BLAZE_ALWAYS_INLINE void zgemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                                const DenseVector<VT2,false>& x, complex<double> alpha, complex<double> beta );

template< typename VT1, typename VT2, typename MT1, bool SO >
BLAZE_ALWAYS_INLINE void zgemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                                const DenseMatrix<MT1,SO>& A, complex<double> alpha, complex<double> beta );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void strmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void strmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void dtrmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void dtrmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void ctrmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void ctrmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void ztrmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
BLAZE_ALWAYS_INLINE void ztrmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
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
// \param y The target left-hand side dense vector.
// \param A The left-hand side dense matrix operand.
// \param x The right-hand side dense vector operand.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for single precision
// operands based on the BLAS cblas_sgemv() function. Note that the function only works for
// vectors and matrices with \c float element type. The attempt to call the function with
// vectors and matrices of any other element type results in a compile time error.
*/
template< typename VT1    // Type of the left-hand side target vector
        , typename MT1    // Type of the left-hand side matrix operand
        , bool SO         // Storage order of the left-hand side matrix operand
        , typename VT2 >  // Type of the right-hand side vector operand
BLAZE_ALWAYS_INLINE void sgemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                                const DenseVector<VT2,false>& x, float alpha, float beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );

   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT2::ElementType );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_sgemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasNoTrans, M, N, alpha,
                (~A).data(), lda, (~x).data(), 1, beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/dense matrix multiplication for single
//        precision operands (\f$ \vec{y}^T=\alpha*\vec{x}^T*A+\beta*\vec{y}^T \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param x The left-hand side dense vector operand.
// \param A The right-hand side dense matrix operand.
// \param alpha The scaling factor for \f$ \vec{x}^T*A \f$.
// \param beta The scaling factor for \f$ \vec{y}^T \f$.
// \return void
//
// This function performs the transpose dense vector/dense matrix multiplication for single
// precision operands based on the BLAS cblas_sgemv() function. Note that the function only
// works for vectors and matrices with \c float element type. The attempt to call the function
// with vectors and matrices of any other element type results in a compile time error.
*/
template< typename VT1  // Type of the left-hand side target vector
        , typename VT2  // Type of the left-hand side vector operand
        , typename MT1  // Type of the right-hand side matrix operand
        , bool SO >     // Storage order of the right-hand side matrix operand
BLAZE_ALWAYS_INLINE void sgemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                                const DenseMatrix<MT1,SO>& A, float alpha, float beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_sgemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasTrans, M, N, alpha,
                (~A).data(), lda, (~x).data(), 1, beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for double precision operands
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
// This function performs the dense matrix/dense vector multiplication for double precision
// operands based on the BLAS cblas_dgemv() function. Note that the function only works for
// vectors and matrices with \c double element type. The attempt to call the function with
// vectors and matrices of any other element type results in a compile time error.
*/
template< typename VT1    // Type of the left-hand side target vector
        , typename MT1    // Type of the left-hand side matrix operand
        , bool SO         // Storage order of the left-hand side matrix operand
        , typename VT2 >  // Type of the right-hand side vector operand
BLAZE_ALWAYS_INLINE void dgemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                                const DenseVector<VT2,false>& x, double alpha, double beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );

   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT2::ElementType );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_dgemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasNoTrans, M, N, alpha,
                (~A).data(), lda, (~x).data(), 1, beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector-dense matrix multiplication for double
//        precision operands (\f$ \vec{y}^T=\alpha*\vec{x}^T*A+\beta*\vec{y}^T \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param x The left-hand side dense vector operand.
// \param A The right-hand side dense matrix operand.
// \param alpha The scaling factor for \f$ \vec{x}^T*A \f$.
// \param beta The scaling factor for \f$ \vec{y}^T \f$.
// \return void
//
// This function performs the transpose dense vector-dense matrix multiplication for double
// precision operands based on the BLAS cblas_dgemv() function. Note that the function only
// works for vectors and matrices with \c double element type. The attempt to call the function
// with vectors and matrices of any other element type results in a compile time error.
*/
template< typename VT1  // Type of the left-hand side target vector
        , typename VT2  // Type of the left-hand side vector operand
        , typename MT1  // Type of the right-hand side matrix operand
        , bool SO >     // Storage order of the right-hand side matrix operand
BLAZE_ALWAYS_INLINE void dgemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                                const DenseMatrix<MT1,SO>& A, double alpha, double beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_dgemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasTrans, M, N, alpha,
                (~A).data(), lda, (~x).data(), 1, beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for single precision complex
//        operands (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The left-hand side dense matrix operand.
// \param x The right-hand side dense vector operand.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for single precision
// complex operands based on the BLAS cblas_cgemv() function. Note that the function only
// works for vectors and matrices with \c complex<float> element type. The attempt to call
// the function with vectors and matrices of any other element type results in a compile
// time error.
*/
template< typename VT1    // Type of the left-hand side target vector
        , typename MT1    // Type of the left-hand side matrix operand
        , bool SO         // Storage order of the left-hand side matrix operand
        , typename VT2 >  // Type of the right-hand side vector operand
BLAZE_ALWAYS_INLINE void cgemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                                const DenseVector<VT2,false>& x, complex<float> alpha, complex<float> beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT2::ElementType::value_type );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_cgemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasNoTrans, M, N, &alpha,
                (~A).data(), lda, (~x).data(), 1, &beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector-dense matrix multiplication for single
//        precision complex operands (\f$ \vec{y}^T=\alpha*\vec{x}^T*A+\beta*\vec{y}^T \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param x The left-hand side dense vector operand.
// \param A The right-hand side dense matrix operand.
// \param alpha The scaling factor for \f$ \vec{x}^T*A \f$.
// \param beta The scaling factor for \f$ \vec{y}^T \f$.
// \return void
//
// This function performs the transpose dense vector-dense matrix multiplication for single
// precision complex operands based on the BLAS cblas_cgemv() function. Note that the function
// only works for vectors and matrices with \c complex<float> element type. The attempt to call
// the function with vectors and matrices of any other element type results in a compile time
// error.
*/
template< typename VT1  // Type of the left-hand side target vector
        , typename VT2  // Type of the left-hand side vector operand
        , typename MT1  // Type of the right-hand side matrix operand
        , bool SO >     // Storage order of the right-hand side matrix operand
BLAZE_ALWAYS_INLINE void cgemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                                const DenseMatrix<MT1,SO>& A, complex<float> alpha, complex<float> beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT2::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_cgemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasTrans, M, N, &alpha,
                (~A).data(), lda, (~x).data(), 1, &beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication for double precision
//        complex operands (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The left-hand side dense matrix operand.
// \param x The right-hand side dense vector operand.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \return void
//
// This function performs the dense matrix/dense vector multiplication for double precision
// complex operands based on the BLAS cblas_zgemv() function. Note that the function only
// works for vectors and matrices with \c complex<double> element type. The attempt to call
// the function with vectors and matrices of any other element type results in a compile
// time error.
*/
template< typename VT1    // Type of the left-hand side target vector
        , typename MT1    // Type of the left-hand side matrix operand
        , bool SO         // Storage order of the left-hand side matrix operand
        , typename VT2 >  // Type of the right-hand side vector operand
BLAZE_ALWAYS_INLINE void zgemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                                const DenseVector<VT2,false>& x, complex<double> alpha, complex<double> beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT2::ElementType::value_type );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_zgemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasNoTrans, M, N, &alpha,
                (~A).data(), lda, (~x).data(), 1, &beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector-dense matrix multiplication for double
//        precision complex operands (\f$ \vec{y}^T=\alpha*\vec{x}^T*A+\beta*\vec{y}^T \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param x The left-hand side dense vector operand.
// \param A The right-hand side dense matrix operand.
// \param alpha The scaling factor for \f$ \vec{x}^T*A \f$.
// \param beta The scaling factor for \f$ \vec{y}^T \f$.
// \return void
//
// This function performs the transpose dense vector-dense matrix multiplication for double
// precision complex operands based on the BLAS cblas_zgemv() function. Note that the function
// only works for vectors and matrices with \c complex<double> element type. The attempt to call
// the function with vectors and matrices of any other element type results in a compile time
// error.
*/
template< typename VT1  // Type of the left-hand side target vector
        , typename VT2  // Type of the left-hand side vector operand
        , typename MT1  // Type of the right-hand side matrix operand
        , bool SO >     // Storage order of the right-hand side matrix operand
BLAZE_ALWAYS_INLINE void zgemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                                const DenseMatrix<MT1,SO>& A, complex<double> alpha, complex<double> beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT2::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_zgemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasTrans, M, N, &alpha,
                (~A).data(), lda, (~x).data(), 1, &beta, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication for single
//        precision operands (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a single precision triangular matrix by a vector
// based on the cblas_strmv() function. Note that the function only works for vectors and matrices
// with \c float element type. The attempt to call the function with vectors and matrices of any
// other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void strmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT::ElementType );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_strmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
                uplo, CblasNoTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/triangular dense matrix multiplication for
//        single precision operands (\f$ \vec{y}^T=\vec{y}^T*A \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a single precision triangular matrix by a vector
// based on the cblas_strmv() function. Note that the function only works for vectors and matrices
// with \c float element type. The attempt to call the function with vectors and matrices of any
// other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void strmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT::ElementType );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_strmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
                uplo, CblasTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication for double
//        precision operands (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a double precision triangular matrix by a vector
// based on the cblas_dtrmv() function. Note that the function only works for vectors and matrices
// with \c double element type. The attempt to call the function with vectors and matrices of any
// other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void dtrmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT::ElementType );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_dtrmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
                uplo, CblasNoTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/triangular dense matrix multiplication for
//        double precision operands (\f$ \vec{y}^T=\vec{y}^T*A \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a double precision triangular matrix by a vector
// based on the cblas_dtrmv() function. Note that the function only works for vectors and matrices
// with \c double element type. The attempt to call the function with vectors and matrices of any
// other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void dtrmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT::ElementType );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_dtrmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
                uplo, CblasTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication for single
//        precision complex operands (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a single precision complex triangular matrix by a
// vector based on the cblas_ctrmv() function. Note that the function only works for vectors and
// matrices with \c complex<float> element type. The attempt to call the function with vectors
// and matrices of any other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void ctrmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT::ElementType::value_type );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_ctrmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
                uplo, CblasNoTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/triangular dense matrix multiplication for
//        single precision complex operands (\f$ \vec{y}^T=\vec{y}^T*A \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a single precision complex triangular matrix by a
// vector based on the cblas_ctrmv() function. Note that the function only works for vectors and
// matrices with \c complex<float> element type. The attempt to call the function with vectors
// and matrices of any other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void ctrmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename VT::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT::ElementType::value_type );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_ctrmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
                uplo, CblasTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication for double
//        precision complex operands (\f$ \vec{y}=A*\vec{y} \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a double precision complex triangular matrix by a
// vector based on the cblas_ztrmv() function. Note that the function only works for vectors and
// matrices with \c complex<double> element type. The attempt to call the function with vectors
// and matrices of any other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void ztrmv( DenseVector<VT,false>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT::ElementType::value_type );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_ztrmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
                uplo, CblasNoTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/triangular dense matrix multiplication for
//        double precision complex operands (\f$ \vec{y}^T=\vec{y}^T*A \f$).
// \ingroup math
//
// \param y The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a double precision complex triangular matrix by a
// vector based on the cblas_ztrmv() function. Note that the function only works for vectors and
// matrices with \c complex<double> element type. The attempt to call the function with vectors
// and matrices of any other element type results in a compile time error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
BLAZE_ALWAYS_INLINE void ztrmv( DenseVector<VT,true>& y, const DenseMatrix<MT,SO>& A,
                                CBLAS_UPLO uplo )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename VT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename VT::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT::ElementType::value_type );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int N  ( numeric_cast<int>( (~A).rows() )    );
   const int lda( numeric_cast<int>( (~A).spacing() ) );

   cblas_ztrmv( ( IsRowMajorMatrix<MT>::value )?( CblasRowMajor ):( CblasColMajor ),
                uplo, CblasTrans, CblasNonUnit, N, (~A).data(), lda, (~y).data(), 1 );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
