//=================================================================================================
/*!
//  \file blaze/math/lapack/gesv.h
//  \brief Header file for LAPACK linear system solver functions (gesv)
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

#ifndef _BLAZE_MATH_LAPACK_GESV_H_
#define _BLAZE_MATH_LAPACK_GESV_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BlasCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/TransposeFlag.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/Exception.h>
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
extern "C" {

void sgesv_( int* n, int* nrhs, float*  A, int* lda, int* ipiv, float*  b, int* ldb, int* info );
void dgesv_( int* n, int* nrhs, double* A, int* lda, int* ipiv, double* b, int* ldb, int* info );
void cgesv_( int* n, int* nrhs, float*  A, int* lda, int* ipiv, float*  b, int* ldb, int* info );
void zgesv_( int* n, int* nrhs, double* A, int* lda, int* ipiv, double* b, int* ldb, int* info );

}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK WRAPPER FUNCTIONS (GESV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK wrapper functions (gesv) */
//@{
void gesv( int* n, int* nrhs, float* A, int* lda, int* ipiv, float* B, int* ldb, int* info );

void gesv( int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info );

void gesv( int* n, int* nrhs, complex<float>* A, int* lda, int* ipiv, complex<float>* B, int* ldb, int* info );

void gesv( int* n, int* nrhs, complex<double>* A, int* lda, int* ipiv, complex<double>* B, int* ldb, int* info );

template< typename MT, typename VT >
void gesv( DenseMatrix<MT,columnMajor>& A, DenseVector<VT,columnVector>& b, int* ipiv );

template< typename MT1, typename MT2 >
void gesv( DenseMatrix<MT1,columnMajor>& A, DenseMatrix<MT2,columnMajor>& B, int* ipiv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a single precision linear system of equations (\f$ A*X=B \f$).
// \ingroup lapack
//
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param B Pointer to the first element of the matrix.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK sgesv() function to compute the solution to the system of linear
// equations \f$ A*X=B \f$, where \a A is a n-by-n matrix and \a X and \a B are n-by-nrhs matrices.
//
// The LU decomposition with partial pivoting and row interchanges is used to factor \a A as

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: \c L is stored in the
// lower part of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L
// are not stored. The factored form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but since factor U(i,i) is exactly
//          singular and the solution could not be computed.
//
// For more information on the sgesv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
void gesv( int* n, int* nrhs, float* A, int* lda, int* ipiv, float* B, int* ldb, int* info )
{
   sgesv_( n, nrhs, A, lda, ipiv, B, ldb, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a double precision linear system of equations (\f$ A*X=B \f$).
// \ingroup lapack
//
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param B Pointer to the first element of the matrix.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK dgesv() function to compute the solution to the system of linear
// equations \f$ A*X=B \f$, where \a A is a n-by-n matrix and \a X and \a B are n-by-nrhs matrices.
//
// The LU decomposition with partial pivoting and row interchanges is used to factor \a A as

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: \c L is stored in the
// lower part of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L
// are not stored. The factored form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but since factor U(i,i) is exactly
//          singular and the solution could not be computed.
//
// For more information on the dgesv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
void gesv( int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info )
{
   dgesv_( n, nrhs, A, lda, ipiv, B, ldb, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a single precision complex linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack
//
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param B Pointer to the first element of the matrix.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK cgesv() function to compute the solution to the system of linear
// equations \f$ A*X=B \f$, where \a A is a n-by-n matrix and \a X and \a B are n-by-nrhs matrices.
//
// The LU decomposition with partial pivoting and row interchanges is used to factor \a A as

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: \c L is stored in the
// lower part of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L
// are not stored. The factored form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but since factor U(i,i) is exactly
//          singular and the solution could not be computed.
//
// For more information on the cgesv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
void gesv( int* n, int* nrhs, complex<float>* A, int* lda, int* ipiv, complex<float>* B, int* ldb, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cgesv_( n, nrhs, reinterpret_cast<float*>( A ), lda, ipiv,
           reinterpret_cast<float*>( B ), ldb, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a double precision complex linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack
//
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param B Pointer to the first element of the matrix.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK zgesv() function to compute the solution to the system of linear
// equations \f$ A*X=B \f$, where \a A is a n-by-n matrix and \a X and \a B are n-by-nrhs matrices.
//
// The LU decomposition with partial pivoting and row interchanges is used to factor \a A as

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: \c L is stored in the
// lower part of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L
// are not stored. The factored form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but since factor U(i,i) is exactly
//          singular and the solution could not be computed.
//
// For more information on the zgesv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
void gesv( int* n, int* nrhs, complex<double>* A, int* lda, int* ipiv, complex<double>* B, int* ldb, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   zgesv_( n, nrhs, reinterpret_cast<double*>( A ), lda, ipiv,
           reinterpret_cast<double*>( B ), ldb, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a linear system of equations (\f$ A*x=b \f$).
// \ingroup lapack
//
// \param A The column-major system matrix.
// \param b The right-hand side vector.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function uses the LAPACK gesv() functions to compute the solution to the system of linear
// equations \f$ A*x=b \f$, where \a A is a column-major n-by-n matrix and \a x and \a b are
// n-dimensional column vectors.
//
// If the function exits successfully, the vector \a b contains the solution of the linear system
// of equations and \a A has been decomposed by means of an LU decomposition with partial pivoting
// and row interchanges. The decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. \c L is stored in the lower part of \a A and \c U is stored in the upper
// part. The unit diagonal elements of \c L are not stored. The factored form of \a A is then
// used to solve the system of equations.
//
// The function fails if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given system matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the gesv() functions (i.e. sgesv(), dgesv(), cgesv(), and zgesv()),
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT, typename VT >
void gesv( DenseMatrix<MT,columnMajor>& A, DenseVector<VT,columnVector>& b, int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( typename MT::ElementType, typename VT::ElementType );

   if( !isSquare( ~A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   int n   ( numeric_cast<int>( (~A).rows() ) );
   int nrhs( 1 );
   int lda ( numeric_cast<int>( (~A).spacing() ) );
   int ldb ( numeric_cast<int>( (~b).size() ) );
   int info( 0 );

   gesv( &n, &nrhs, (~A).data(), &lda, ipiv, (~b).data(), &ldb, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid function argument" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a linear system of equations (\f$ A*X=B \f$).
// \ingroup lapack
//
// \param A The system matrix.
// \param B The matrix of right-hand sides.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function uses the LAPACK gesv() functions to compute the solution to the system of linear
// equations \f$ A*X=B \f$, where \a A is a column-major n-by-n matrix and \a X and \a B are
// column-major n-by-m matrices.
//
// If the function exits successfully, the matrix \a B contains the solutions of the linear system
// of equations and \a A has been decomposed by means of an LU decomposition with partial pivoting
// and row interchanges. The decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. \c L is stored in the lower part of \a A and \c U is stored in the upper
// part. The unit diagonal elements of \c L are not stored. The factored form of \a A is then
// used to solve the system of equations.
//
// The function fails if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given system matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the gesv() functions (i.e. sgesv(), dgesv(), cgesv(), and zgesv()),
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT1, typename MT2 >
void gesv( DenseMatrix<MT1,columnMajor>& A, DenseMatrix<MT2,columnMajor>& B, int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT1::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( typename MT1::ElementType, typename MT2::ElementType );

   if( !isSquare( ~A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   int n   ( numeric_cast<int>( (~A).rows()    ) );
   int nrhs( numeric_cast<int>( (~B).columns() ) );
   int lda ( numeric_cast<int>( (~A).spacing() ) );
   int ldb ( numeric_cast<int>( (~B).spacing() ) );
   int info( 0 );

   gesv( &n, &nrhs, (~A).data(), &lda, ipiv, (~B).data(), &ldb, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid function argument" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
