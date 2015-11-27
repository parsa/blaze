//=================================================================================================
/*!
//  \file blaze/math/lapack/Inversion.h
//  \brief Header file for LAPACK matrix inversion functions
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

#ifndef _BLAZE_MATH_LAPACK_INVERSION_H_
#define _BLAZE_MATH_LAPACK_INVERSION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/lapack/LU.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/Complex.h>
#include <blaze/util/constraints/Double.h>
#include <blaze/util/constraints/Float.h>
#include <blaze/util/Exception.h>
#include <blaze/util/UniqueArray.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

extern "C" {

void sgetri_( int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info );
void dgetri_( int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info );
void cgetri_( int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info );
void zgetri_( int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info );

}




//=================================================================================================
//
//  LAPACK INVERSION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK inversion functions */
//@{
template< typename MT, bool SO >
inline void sgetri( DenseMatrix<MT,SO>& A, const int* ipiv );

template< typename MT, bool SO >
inline void dgetri( DenseMatrix<MT,SO>& A, const int* ipiv );

template< typename MT, bool SO >
inline void cgetri( DenseMatrix<MT,SO>& A, const int* ipiv );

template< typename MT, bool SO >
inline void zgetri( DenseMatrix<MT,SO>& A, const int* ipiv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense single precision matrix.
// \ingroup lapack
//
// \param A The matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function performs the dense matrix inversion based on the LAPACK sgetri() function for
// single precision matrices that have already been factorized by sgetrf() function. Note that the
// function only works for general, non-adapted matrices with \c float element type. The attempt
// to call the function with adaptors or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the sgetri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown, \c A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void sgetri( DenseMatrix<MT,SO>& A, const int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT::ElementType );

   if( !IsSquare<MT>::value && !isSquare( A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   int n    ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda  ( boost::numeric_cast<int>( (~A).spacing() ) );
   int lwork( n*lda );
   int info ( 0 );

   const UniqueArray<float> work( new float[lwork] );

   sgetri_( &n, (~A).data(), &lda, const_cast<int*>( ipiv ), work.get(), &lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for matrix inversion" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense double precision matrix.
// \ingroup lapack
//
// \param A The matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function performs the dense matrix inversion based on the LAPACK dgetri() function for
// double precision matrices that have already been factorized by dgetrf() function. Note that the
// function only works for general, non-adapted matrices with \c double element type. The attempt
// to call the function with adaptors or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the dgetri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown, \c A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void dgetri( DenseMatrix<MT,SO>& A, const int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT::ElementType );

   if( !IsSquare<MT>::value && !isSquare( A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   int n    ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda  ( boost::numeric_cast<int>( (~A).spacing() ) );
   int lwork( n*lda );
   int info ( 0 );

   const UniqueArray<double> work( new double[lwork] );

   dgetri_( &n, (~A).data(), &lda, const_cast<int*>( ipiv ), work.get(), &lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for matrix inversion" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense single precision complex matrix.
// \ingroup lapack
//
// \param A The matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function performs the dense matrix inversion based on the LAPACK cgetri() function for
// single precision complex matrices that have already been factorized by cgetrf() function. Note
// that the function only works for general, non-adapted matrices with \c complex<float> element
// type. The attempt to call the function with adaptors or matrices of any other element type
// results in a compile time error!
//
// The function fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the cgetri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown, \c A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void cgetri( DenseMatrix<MT,SO>& A, const int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE  ( typename MT::ElementType::value_type );

   if( !IsSquare<MT>::value && !isSquare( A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   int n    ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda  ( boost::numeric_cast<int>( (~A).spacing() ) );
   int lwork( n*lda );
   int info ( 0 );

   const UniqueArray< complex<float> > work( new complex<float>[lwork] );

   cgetri_( &n, reinterpret_cast<float*>( (~A).data() ), &lda, const_cast<int*>( ipiv ),
            reinterpret_cast<float*>( work.get() ), &lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for matrix inversion" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense double precision complex matrix.
// \ingroup lapack
//
// \param A The matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function performs the dense matrix inversion based on the LAPACK zgetri() function for
// double precision complex matrices that have already been factorized by zgetrf() function. Note
// that the function only works for general, non-adapted matrices with \c complex<double> element
// type. The attempt to call the function with adaptors or matrices of any other element type
// results in a compile time error!
//
// The function fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// For more information on the zgetri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown, \c A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void zgetri( DenseMatrix<MT,SO>& A, const int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE ( typename MT::ElementType::value_type );

   if( !IsSquare<MT>::value && !isSquare( A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   int n    ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda  ( boost::numeric_cast<int>( (~A).spacing() ) );
   int lwork( n*lda );
   int info ( 0 );

   const UniqueArray< complex<double> > work( new complex<double>[lwork] );

   zgetri_( &n, reinterpret_cast<double*>( (~A).data() ), &lda, const_cast<int*>( ipiv ),
            reinterpret_cast<double*>( work.get() ), &lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for matrix inversion" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
