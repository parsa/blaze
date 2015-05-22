//=================================================================================================
/*!
//  \file blaze/math/blas/Level3.h
//  \brief Header file for BLAS level 3 functions
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

#ifndef _BLAZE_MATH_BLAS_LEVEL3_H_
#define _BLAZE_MATH_BLAS_LEVEL3_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/ConstDataAccess.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSymmetric.h>
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
//  BLAS LEVEL 3 FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS level 3 functions */
//@{
#if BLAZE_BLAS_MODE

template< typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3 >
BLAZE_ALWAYS_INLINE void sgemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                                const DenseMatrix<MT3,SO3>& B, float alpha, float beta );

template< typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3 >
BLAZE_ALWAYS_INLINE void dgemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                                const DenseMatrix<MT3,SO3>& B, double alpha, double beta );

template< typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3 >
BLAZE_ALWAYS_INLINE void cgemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                                const DenseMatrix<MT3,SO3>& B, complex<float> alpha, complex<float> beta );

template< typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3 >
BLAZE_ALWAYS_INLINE void zgemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                                const DenseMatrix<MT3,SO3>& B, complex<double> alpha, complex<double> beta );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
BLAZE_ALWAYS_INLINE void strmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                                CBLAS_SIDE side, CBLAS_UPLO uplo, float alpha );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
BLAZE_ALWAYS_INLINE void dtrmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                                CBLAS_SIDE side, CBLAS_UPLO uplo, double alpha );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
BLAZE_ALWAYS_INLINE void ctrmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                                CBLAS_SIDE side, CBLAS_UPLO uplo, complex<float> alpha );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
BLAZE_ALWAYS_INLINE void ztrmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                                CBLAS_SIDE side, CBLAS_UPLO uplo, complex<double> alpha );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication with single precision
//        matrices (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup math
//
// \param C The target left-hand side dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication for single precision
// matrices based on the BLAS cblas_sgemm() and cblas_ssymm() functions. Note that the function
// only works for matrices with \c float element type. The attempt to call the function with
// matrices of any other element type results in a compile time error.
*/
template< typename MT1  // Type of the left-hand side target matrix
        , bool SO1      // Storage order of the left-hand side target matrix
        , typename MT2  // Type of the left-hand side matrix operand
        , bool SO2      // Storage order of the left-hand side matrix operand
        , typename MT3  // Type of the right-hand side matrix operand
        , bool SO3 >    // Storage order of the right-hand side matrix operand
BLAZE_ALWAYS_INLINE void sgemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                                const DenseMatrix<MT3,SO3>& B, float alpha, float beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT3::ElementType );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~B).columns() ) );
   const int K  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );
   const int ldb( numeric_cast<int>( (~B).spacing() ) );
   const int ldc( numeric_cast<int>( (~C).spacing() ) );

   if( IsSymmetric<MT2>::value && ( SO1 == SO3 ) ) {
      cblas_ssymm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   CblasLeft,
                   ( IsRowMajorMatrix<MT2>::value )?( CblasLower ):( CblasUpper ),
                   M, N, alpha, (~A).data(), lda, (~B).data(), ldb, beta, (~C).data(), ldc );
   }
   else if( IsSymmetric<MT3>::value && ( SO1 == SO2 ) ) {
      cblas_ssymm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   CblasRight,
                   ( IsRowMajorMatrix<MT3>::value )?( CblasLower ):( CblasUpper ),
                   M, N, alpha, (~B).data(), ldb, (~A).data(), lda, beta, (~C).data(), ldc );
   }
   else {
      cblas_sgemm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
                   ( SO1 == SO3 )?( CblasNoTrans ):( CblasTrans ),
                   M, N, K, alpha, (~A).data(), lda, (~B).data(), ldb, beta, (~C).data(), ldc );
   }
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication with double precision
//        matrices (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup math
//
// \param C The target left-hand side dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication for double precision
// matrices based on the BLAS cblas_dgemm() and cblas_dsymm() functions. Note that the function
// only works for matrices with \c double element type. The attempt to call the function with
// matrices of any other element type results in a compile time error.
*/
template< typename MT1  // Type of the left-hand side target matrix
        , bool SO1      // Storage order of the left-hand side target matrix
        , typename MT2  // Type of the left-hand side matrix operand
        , bool SO2      // Storage order of the left-hand side matrix operand
        , typename MT3  // Type of the right-hand side matrix operand
        , bool SO3 >    // Storage order of the right-hand side matrix operand
BLAZE_ALWAYS_INLINE void dgemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                                const DenseMatrix<MT3,SO3>& B, double alpha, double beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT3::ElementType );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~B).columns() ) );
   const int K  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );
   const int ldb( numeric_cast<int>( (~B).spacing() ) );
   const int ldc( numeric_cast<int>( (~C).spacing() ) );

   if( IsSymmetric<MT2>::value && ( SO1 == SO3 ) ) {
      cblas_dsymm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   CblasLeft,
                   ( IsRowMajorMatrix<MT2>::value )?( CblasLower ):( CblasUpper ),
                   M, N, alpha, (~A).data(), lda, (~B).data(), ldb, beta, (~C).data(), ldc );
   }
   else if( IsSymmetric<MT3>::value && ( SO1 == SO2 ) ) {
      cblas_dsymm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   CblasRight,
                   ( IsRowMajorMatrix<MT3>::value )?( CblasLower ):( CblasUpper ),
                   M, N, alpha, (~B).data(), ldb, (~A).data(), lda, beta, (~C).data(), ldc );
   }
   else {
      cblas_dgemm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
                   ( SO1 == SO3 )?( CblasNoTrans ):( CblasTrans ),
                   M, N, K, alpha, (~A).data(), lda, (~B).data(), ldb, beta, (~C).data(), ldc );
   }
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication with single precision
//        complex matrices (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup math
//
// \param C The target left-hand side dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication for single precision
// complex matrices based on the BLAS cblas_cgemm() and cblas_csymm() functions. Note that
// the function only works for matrices with \c complex<float> element type. The attempt to
// call the function with matrices of any other element type results in a compile time error.
*/
template< typename MT1  // Type of the left-hand side target matrix
        , bool SO1      // Storage order of the left-hand side target matrix
        , typename MT2  // Type of the left-hand side matrix operand
        , bool SO2      // Storage order of the left-hand side matrix operand
        , typename MT3  // Type of the right-hand side matrix operand
        , bool SO3 >    // Storage order of the right-hand side matrix operand
BLAZE_ALWAYS_INLINE void cgemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                                const DenseMatrix<MT3,SO3>& B, complex<float> alpha, complex<float> beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT3::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT2::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT3::ElementType::value_type );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~B).columns() ) );
   const int K  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );
   const int ldb( numeric_cast<int>( (~B).spacing() ) );
   const int ldc( numeric_cast<int>( (~C).spacing() ) );

   if( IsSymmetric<MT2>::value && ( SO1 == SO3 ) ) {
      cblas_csymm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   CblasLeft,
                   ( IsRowMajorMatrix<MT2>::value )?( CblasLower ):( CblasUpper ),
                   M, N, &alpha, (~A).data(), lda, (~B).data(), ldb, &beta, (~C).data(), ldc );
   }
   else if( IsSymmetric<MT3>::value && ( SO1 == SO2 ) ) {
      cblas_csymm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   CblasRight,
                   ( IsRowMajorMatrix<MT3>::value )?( CblasLower ):( CblasUpper ),
                   M, N, &alpha, (~B).data(), ldb, (~A).data(), lda, &beta, (~C).data(), ldc );
   }
   else {
      cblas_cgemm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
                   ( SO1 == SO3 )?( CblasNoTrans ):( CblasTrans ),
                   M, N, K, &alpha, (~A).data(), lda, (~B).data(), ldb, &beta, (~C).data(), ldc );
   }
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication with double precision
//        complex matrices (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup math
//
// \param C The target left-hand side dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication for double precision
// complex matrices based on the BLAS cblas_zgemm() and cblas_zsymm() functions. Note that
// the function only works for matrices with \c complex<double> element type. The attempt to
// call the function with matrices of any other element type results in a compile time error.
*/
template< typename MT1  // Type of the left-hand side target matrix
        , bool SO1      // Storage order of the left-hand side target matrix
        , typename MT2  // Type of the left-hand side matrix operand
        , bool SO2      // Storage order of the left-hand side matrix operand
        , typename MT3  // Type of the right-hand side matrix operand
        , bool SO3 >    // Storage order of the right-hand side matrix operand
BLAZE_ALWAYS_INLINE void zgemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                                const DenseMatrix<MT3,SO3>& B, complex<double> alpha, complex<double> beta )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT3::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT2::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT3::ElementType::value_type );

   const int M  ( numeric_cast<int>( (~A).rows() )    );
   const int N  ( numeric_cast<int>( (~B).columns() ) );
   const int K  ( numeric_cast<int>( (~A).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );
   const int ldb( numeric_cast<int>( (~B).spacing() ) );
   const int ldc( numeric_cast<int>( (~C).spacing() ) );

   if( IsSymmetric<MT2>::value && ( SO1 == SO3 ) ) {
      cblas_zsymm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   CblasLeft,
                   ( IsRowMajorMatrix<MT2>::value )?( CblasLower ):( CblasUpper ),
                   M, N, &alpha, (~A).data(), lda, (~B).data(), ldb, &beta, (~C).data(), ldc );
   }
   else if( IsSymmetric<MT3>::value && ( SO1 == SO2 ) ) {
      cblas_zsymm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   CblasRight,
                   ( IsRowMajorMatrix<MT3>::value )?( CblasLower ):( CblasUpper ),
                   M, N, &alpha, (~B).data(), ldb, (~A).data(), lda, &beta, (~C).data(), ldc );
   }
   else {
      cblas_zgemm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                   ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
                   ( SO1 == SO3 )?( CblasNoTrans ):( CblasTrans ),
                   M, N, K, &alpha, (~A).data(), lda, (~B).data(), ldb, &beta, (~C).data(), ldc );
   }
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense matrix multiplication with single
//        precision matrices (\f$ B=\alpha*A*B \f$ or \f$ B=\alpha*B*A \f$).
// \ingroup math
//
// \param B The target dense matrix.
// \param A The dense matrix multiplication operand.
// \param side \a CblasLeft to compute \f$ B=\alpha*A*B \f$, \a CblasRight to compute \f$ B=\alpha*B*A \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param alpha The scaling factor for \f$ A*B \f$ or \f$ B*A \f$.
// \return void
//
// This function performs the scaling and multiplication of a triangular matrix by a matrix
// based on the cblas_strmm() function. Note that the function only works for matrices with
// \c float element type. The attempt to call the function with matrices of any other element
// type results in a compile time error. Also, matrix \a A is expected to be a square matrix.
*/
template< typename MT1  // Type of the left-hand side target matrix
        , bool SO1      // Storage order of the left-hand side target matrix
        , typename MT2  // Type of the left-hand side matrix operand
        , bool SO2 >    // Storage order of the left-hand side matrix operand
BLAZE_ALWAYS_INLINE void strmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                                CBLAS_SIDE side, CBLAS_UPLO uplo, float alpha )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT2::ElementType );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( side == CblasLeft  || side == CblasRight, "Invalid side argument detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int M  ( numeric_cast<int>( (~B).rows() )    );
   const int N  ( numeric_cast<int>( (~B).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );
   const int ldb( numeric_cast<int>( (~B).spacing() ) );

   cblas_strmm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                side,
                ( SO1 == SO2 )?( uplo )
                              :( ( uplo == CblasLower )?( CblasUpper ):( CblasLower ) ),
                ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
                CblasNonUnit,
                M, N, alpha, (~A).data(), lda, (~B).data(), ldb );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense matrix multiplication with double
//        precision matrices (\f$ B=\alpha*A*B \f$ or \f$ B=\alpha*B*A \f$).
// \ingroup math
//
// \param B The target dense matrix.
// \param A The dense matrix multiplication operand.
// \param side \a CblasLeft to compute \f$ B=\alpha*A*B \f$, \a CblasRight to compute \f$ B=\alpha*B*A \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param alpha The scaling factor for \f$ A*B \f$ or \f$ B*A \f$.
// \return void
//
// This function performs the scaling and multiplication of a triangular matrix by a matrix
// based on the cblas_dtrmm() function. Note that the function only works for matrices with
// \c double element type. The attempt to call the function with matrices of any other element
// type results in a compile time error. Also, matrix \a A is expected to be a square matrix.
*/
template< typename MT1  // Type of the left-hand side target matrix
        , bool SO1      // Storage order of the left-hand side target matrix
        , typename MT2  // Type of the left-hand side matrix operand
        , bool SO2 >    // Storage order of the left-hand side matrix operand
BLAZE_ALWAYS_INLINE void dtrmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                                CBLAS_SIDE side, CBLAS_UPLO uplo, double alpha )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT2::ElementType );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( side == CblasLeft  || side == CblasRight, "Invalid side argument detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int M  ( numeric_cast<int>( (~B).rows() )    );
   const int N  ( numeric_cast<int>( (~B).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );
   const int ldb( numeric_cast<int>( (~B).spacing() ) );

   cblas_dtrmm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                side,
                ( SO1 == SO2 )?( uplo )
                              :( ( uplo == CblasLower )?( CblasUpper ):( CblasLower ) ),
                ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
                CblasNonUnit,
                M, N, alpha, (~A).data(), lda, (~B).data(), ldb );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense matrix multiplication with single
//        precision complex matrices (\f$ B=\alpha*A*B \f$ or \f$ B=\alpha*B*A \f$).
// \ingroup math
//
// \param B The target dense matrix.
// \param A The dense matrix multiplication operand.
// \param side \a CblasLeft to compute \f$ B=\alpha*A*B \f$, \a CblasRight to compute \f$ B=\alpha*B*A \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param alpha The scaling factor for \f$ A*B \f$ or \f$ B*A \f$.
// \return void
//
// This function performs the scaling and multiplication of a triangular matrix by a matrix
// based on the cblas_ctrmm() function. Note that the function only works for matrices with
// \c complex<float> element type. The attempt to call the function with matrices of any
// other element type results in a compile time error. Also, matrix \a A is expected to be
// a square matrix.
*/
template< typename MT1  // Type of the left-hand side target matrix
        , bool SO1      // Storage order of the left-hand side target matrix
        , typename MT2  // Type of the left-hand side matrix operand
        , bool SO2 >    // Storage order of the left-hand side matrix operand
BLAZE_ALWAYS_INLINE void ctrmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                                CBLAS_SIDE side, CBLAS_UPLO uplo, complex<float> alpha )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT2::ElementType::value_type );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( side == CblasLeft  || side == CblasRight, "Invalid side argument detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int M  ( numeric_cast<int>( (~B).rows() )    );
   const int N  ( numeric_cast<int>( (~B).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );
   const int ldb( numeric_cast<int>( (~B).spacing() ) );

   cblas_ctrmm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                side,
                ( SO1 == SO2 )?( uplo )
                              :( ( uplo == CblasLower )?( CblasUpper ):( CblasLower ) ),
                ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
                CblasNonUnit,
                M, N, &alpha, (~A).data(), lda, (~B).data(), ldb );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense matrix multiplication with double
//        precision complex matrices (\f$ B=\alpha*A*B \f$ or \f$ B=\alpha*B*A \f$).
// \ingroup math
//
// \param B The target dense matrix.
// \param A The dense matrix multiplication operand.
// \param side \a CblasLeft to compute \f$ B=\alpha*A*B \f$, \a CblasRight to compute \f$ B=\alpha*B*A \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param alpha The scaling factor for \f$ A*B \f$ or \f$ B*A \f$.
// \return void
//
// This function performs the scaling and multiplication of a triangular matrix by a matrix
// based on the cblas_ztrmm() function. Note that the function only works for matrices with
// \c complex<double> element type. The attempt to call the function with matrices of any
// other element type results in a compile time error. Also, matrix \a A is expected to be
// a square matrix.
*/
template< typename MT1  // Type of the left-hand side target matrix
        , bool SO1      // Storage order of the left-hand side target matrix
        , typename MT2  // Type of the left-hand side matrix operand
        , bool SO2 >    // Storage order of the left-hand side matrix operand
BLAZE_ALWAYS_INLINE void ztrmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                                CBLAS_SIDE side, CBLAS_UPLO uplo, complex<double> alpha )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT2::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT1::ElementType::value_type );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT2::ElementType::value_type );

   BLAZE_INTERNAL_ASSERT( (~A).rows() == (~A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( side == CblasLeft  || side == CblasRight, "Invalid side argument detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const int M  ( numeric_cast<int>( (~B).rows() )    );
   const int N  ( numeric_cast<int>( (~B).columns() ) );
   const int lda( numeric_cast<int>( (~A).spacing() ) );
   const int ldb( numeric_cast<int>( (~B).spacing() ) );

   cblas_ztrmm( ( IsRowMajorMatrix<MT1>::value )?( CblasRowMajor ):( CblasColMajor ),
                side,
                ( SO1 == SO2 )?( uplo )
                              :( ( uplo == CblasLower )?( CblasUpper ):( CblasLower ) ),
                ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
                CblasNonUnit,
                M, N, &alpha, (~A).data(), lda, (~B).data(), ldb );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
