//=================================================================================================
/*!
//  \file blaze/math/lapack/gges.h
//  \brief Header file for the LAPACK general matrix eigenvalue functions (gges)
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

#ifndef _BLAZE_MATH_LAPACK_GGES_H_
#define _BLAZE_MATH_LAPACK_GGES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <memory>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/Contiguous.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/lapack/clapack/gges.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/constraints/Complex.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/NumericCast.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK GENERALIZED MATRIX EIGENVALUE FUNCTIONS (GGES)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK generalized Schur factorization functions (gges) */
//@{
template< typename MT1, bool SO1, typename MT2, bool SO2, 
   typename VT1, bool TF1, typename VT2, bool TF2,
   typename MT3, bool SO3, typename MT4, bool SO4 >
inline void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
   DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta, 
   DenseMatrix<MT3,SO3>& VSL, DenseMatrix<MT4,SO4>& VSR );


template<
   typename MT1, bool SO1, typename MT2, bool SO2, 
   typename VT1, bool TF1, typename VT2, bool TF2,
   typename MT3, bool SO3, typename MT4, bool SO4 >
inline void gges( int (*selctg)(ElementType_t<VT2> const*, ElementType_t<VT2> const*, ElementType_t<VT2> const*), 
   DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
   DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta, 
   DenseMatrix<MT3,SO3>& VSL, DenseMatrix<MT4,SO4>& VSR );

//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur factorization
// of the given pair of dense general matrices.
// \ingroup lapack_eigenvalue
//
// \param A On entry, the first of the pair of matrices. On exit, \a A has been overwritten 
// by its generalized Schur form S.
// \param B On entry, the second of the pair of matrices. On exit, \a B has been overwritten
// by its generalized Schur form T.
// \param alpha The resulting complex vector of eigenvalues numerator. Resized if necessary.
// \param beta The resulting real vector of eigenvalues denominator. Resized if necessary.
// \param VSL The matrix of resulting left Schur vectors. Resized if necessary.
// \param VSR The matrix of resulting right Schur vectors. Resized if necessary.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Vector or matrix cannot be resized.
// \exception std::runtime_error Schur factorization computation failed.
//
// This function computes for a pair of N-by-N real nonsymmetric matrices (A,B),
//  the generalized eigenvalues, the generalized real Schur form (S,T),
//  and the left and right matrices of Schur vectors (VSL and
//  VSR). This gives the generalized Schur factorization
// 
//           (A^FA,B^FB) = ( (VSL^FL)*(S^FA)*(VSR^FR)^T, (VSL^FL)*T*(VSR^FR)^T )
//
// where FA, FB, FL, FR are transposition flags:
//  FA = 1 if \a A is column-major and FA = T (transpose) if \a A is row-major,
//  FB = 1 if \a B is column-major and FB = T (transpose) if \a B is row-major,
//  FL = 1 if \a VSL is column-major and FL = T (transpose) if \a VSL is row-major,
//  FR = 1 if \a VSR is column-major and FR = T (transpose) if \a VSR is row-major.
// 
// A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
//  or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
//  usually represented as the pair (alpha,beta), as there is a
//  reasonable interpretation for beta=0 or both being zero.
//  The complex eigenvalues are returned as numerators and denominators in the given vectors
// \a alpha, \a beta, which are resized to the correct size (if possible and necessary).
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given matrix \a B is not a square matrix;
//  - ... the size of the given matrces \a A and \a B don't match;
//  - ... the given vector \a alpha is a fixed size vector and the size doesn't match;
//  - ... the given vector \a beta is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a VSL is a fixed size matrix and the size doesn't match;
//  - ... the given matrix \a VSR is a fixed size matrix and the size doesn't match;
//  - ... the Schur factorization computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL ), B( 5UL, 5UL );  // The general matrces A, B
   // ... Initialization

   DynamicVector<complex<double>,columnVector> alpha( 5UL );  // The numerator vector for the complex eigenvalues
   DynamicVector<double,columnVector> beta( 5UL );  // The denominator vector for the complex eigenvalues
   DynamicMatrix<double,columnVector> VSL( 5UL, 5UL ), VSR( 5UL, 5UL );  // The matrices of left and right Schur vectors

   gges( A, B, alpha, beta, VSL, VSR );
   \endcode
//
// For more information on the gges() functions (i.e. sgges(), dgges(), cgges(), and zgges())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< 
   typename MT1, bool SO1, typename MT2, bool SO2, 
   typename VT1, bool TF1, typename VT2, bool TF2,
   typename MT3, bool SO3, typename MT4, bool SO4 >
inline void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
   DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta, 
   DenseMatrix<MT3,SO3>& VSL, DenseMatrix<MT4,SO4>& VSR )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT4> );
   
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPLEX_TYPE( ElementType_t<VT2> );

   using RT = ElementType_t<VT2>;

   const size_t N( (~A).rows() );

   if( !isSquare( ~A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( ~B, N, N, false );
   resize( ~alpha, N, false );
   resize( ~beta, N, false );
   resize( ~VSL, N, N, false );
   resize( ~VSR, N, N, false );

   if( N == 0UL ) {
      return;
   }

   gges_backend( (int (*)(RT const*, RT const*, RT const*))nullptr, ~A, ~B, ~alpha, ~beta, ~VSL, ~VSR );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur factorization
// of the given pair of dense general matrices with eigenvalue selection.
// \ingroup lapack_eigenvalue
//
// \param selctg is a function of three real arguments.
// \a selctg is used to select eigenvalues to sort to the top left of the Schur form.
// An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
// selctg(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
// one of a complex conjugate pair of eigenvalues is selected,
// then both complex eigenvalues are selected.
// \param A On entry, the first of the pair of matrices. On exit, \a A has been overwritten 
// by its generalized Schur form S.
// \param B On entry, the second of the pair of matrices. On exit, \a B has been overwritten
// by its generalized Schur form T.
// \param alpha The resulting complex vector of eigenvalues numerator. Resized if necessary.
// \param beta The resulting real vector of eigenvalues denominator. Resized if necessary.
// \param VSL The matrix of resulting left Schur vectors. Resized if necessary.
// \param VSR The matrix of resulting right Schur vectors. Resized if necessary.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Vector or matrix cannot be resized.
// \exception std::runtime_error Schur factorization computation failed.
//
// This function computes for a pair of N-by-N real nonsymmetric matrices (A,B),
//  the generalized eigenvalues, the generalized real Schur form (S,T),
//  and the left and right matrices of Schur vectors (VSL and
//  VSR). This gives the generalized Schur factorization
// 
//           (A^FA,B^FB) = ( (VSL^FL)*(S^FA)*(VSR^FR)^T, (VSL^FL)*T*(VSR^FR)^T )
//
// where FA, FB, FL, FR are transposition flags:
//  FA = 1 if \a A is column-major and FA = T (transpose) if \a A is row-major,
//  FB = 1 if \a B is column-major and FB = T (transpose) if \a B is row-major,
//  FL = 1 if \a VSL is column-major and FL = T (transpose) if \a VSL is row-major,
//  FR = 1 if \a VSR is column-major and FR = T (transpose) if \a VSR is row-major.
// 
// A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
//  or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
//  usually represented as the pair (alpha,beta), as there is a
//  reasonable interpretation for beta=0 or both being zero.
//  The complex eigenvalues are returned as numerators and denominators in the given vectors
// \a alpha, \a beta, which are resized to the correct size (if possible and necessary).
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given matrix \a B is not a square matrix;
//  - ... the size of the given matrces \a A and \a B don't match;
//  - ... the given vector \a alpha is a fixed size vector and the size doesn't match;
//  - ... the given vector \a beta is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a VSL is a fixed size matrix and the size doesn't match;
//  - ... the given matrix \a VSR is a fixed size matrix and the size doesn't match;
//  - ... the Schur factorization computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL ), B( 5UL, 5UL );  // The general matrces A, B
   // ... Initialization

   DynamicVector<complex<double>,columnVector> alpha( 5UL );  // The numerator vector for the complex eigenvalues
   DynamicVector<double,columnVector> beta( 5UL );  // The denominator vector for the complex eigenvalues
   DynamicMatrix<double,columnVector> VSL( 5UL, 5UL ), VSR( 5UL, 5UL );  // The matrices of left and right Schur vectors

   gges( A, B, alpha, beta, VSL, VSR );
   \endcode
//
// For more information on the gges() functions (i.e. sgges(), dgges(), cgges(), and zgges())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< 
   typename MT1, bool SO1, typename MT2, bool SO2, 
   typename VT1, bool TF1, typename VT2, bool TF2,
   typename MT3, bool SO3, typename MT4, bool SO4 >
inline void gges(
   int (*selctg)(ElementType_t<VT2> const*, ElementType_t<VT2> const*, ElementType_t<VT2> const*),
   DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
   DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta, 
   DenseMatrix<MT3,SO3>& VSL, DenseMatrix<MT4,SO4>& VSR )
{
   using RT = ElementType_t<VT2>;
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPLEX_TYPE( RT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT4> );
   
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPLEX_TYPE( ElementType_t<VT2> );

   const size_t N( (~A).rows() );

   if( !isSquare( ~A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( ~B, N, N, false );
   resize( ~alpha, N, false );
   resize( ~beta, N, false );
   resize( ~VSL, N, N, false );
   resize( ~VSR, N, N, false );

   if( N == 0UL ) {
      return;
   }

   gges_backend( selctg, ~A, ~B, ~alpha, ~beta, ~VSL, ~VSR );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gges kernel for real general matrices.
// \ingroup lapack_eigenvalue
//
// This function is the backend implementation for computing the generalized Schur factorization
// of the given pair of real dense general matrices.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gges() function.
*/
template< 
   typename MT1, bool SO1, typename MT2, bool SO2, 
   typename VT1, bool TF1, typename VT2, bool TF2,
   typename MT3, bool SO3, typename MT4, bool SO4 >
inline auto gges_backend( int (*selctg)(ElementType_t<VT2> const*, ElementType_t<VT2> const*, ElementType_t<VT2> const*), 
   DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, 
   DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta, 
   DenseMatrix<MT3,SO3>& VSL, DenseMatrix<MT4,SO4>& VSR )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   using RT = ElementType_t<VT2>;
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPLEX_TYPE( RT );

   int const n = numeric_cast<int>( rows( A ) );
   
   BLAZE_INTERNAL_ASSERT( isSquare( ~A ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( ~B ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( ~VSL ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( ~VSR ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( rows( B ) == n, "Invalid matrix size detected" );
   BLAZE_INTERNAL_ASSERT( rows( VSL ) == n, "Invalid matrix size detected" );
   BLAZE_INTERNAL_ASSERT( rows( VSR ) == n, "Invalid matrix size detected" );
   
   BLAZE_INTERNAL_ASSERT( size( alphar ) == n, "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( size( alphal ) == n, "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( size( beta ) == n, "Invalid vector dimension detected" );

   using CT = ElementType_t<VT1>;
   using BT = ElementType_t<MT1>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   int const lda = numeric_cast<int>( spacing( A ) );
   int const ldb = numeric_cast<int>( spacing( B ) );
   int const ldvsl = numeric_cast<int>( spacing( VSL ) );
   int const ldvsr = numeric_cast<int>( spacing( VSR ) );
   int info = 0;
   int sdim = 0;

   int lwork = std::max( 8*n, 6*n + 16 );
   const std::unique_ptr<BT[]> alphar( new BT[n] );
   const std::unique_ptr<BT[]> alphai( new BT[n] );
   const std::unique_ptr<BT[]> work( new BT[std::max(1, lwork)] );
   const std::unique_ptr<int[]> bwork( new int[n] );

   gges( 'V', 'V', selctg ? 'S' : 'N', (int (*)(RT*, RT*, RT*))selctg, n,
      data( A ), lda, data( B ), ldb, &sdim, 
      alphar.get(), alphai.get(), data( beta ), data( VSL ), ldvsl, data( VSR ), ldvsr, 
      work.get(), lwork, bwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for generalized eigenvalue decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Generalized eigenvalue decomposition failed" );
   }

   for( size_t i=0UL; i < rows( A ); ++i ) {
      (~alpha)[i] = CT( alphar[i], alphai[i] );
   }
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
