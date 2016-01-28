//=================================================================================================
/*!
//  \file blaze/math/dense/QR.h
//  \brief Header file for the dense matrix in-place QR decomposition
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

#ifndef _BLAZE_MATH_DENSE_QR_H_
#define _BLAZE_MATH_DENSE_QR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BlasCompatible.h>
#include <blaze/math/constraints/Hermitian.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/StrictlyTriangular.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/Functions.h>
#include <blaze/math/lapack/geqrf.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/views/Column.h>
#include <blaze/math/views/DenseColumn.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/Exception.h>
#include <blaze/util/mpl/If.h>


namespace blaze {

//=================================================================================================
//
//  QR DECOMPOSITION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name QR decomposition functions */
//@{
template< typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3 >
void qr( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& Q, DenseMatrix<MT3,SO3>& R );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief QR decomposition of the given dense matrix.
// \ingroup dense_matrix
//
// \param A The matrix to be decomposed.
// \param Q The resulting Q matrix.
// \param R The resulting R matrix.
// \return void
// \exception std::invalid_argument Dimensions of fixed size matrix do not match.
//
// This function performs the dense matrix QR decomposition of a general m-by-n matrix. The
// resulting decomposition has the form

                              \f[ A = Q \cdot R, \f]

// where \c Q is a general m-by-m matrix and \c R is an upper trapezoidal m-by-n matrix. The
// decomposition is written to the two distinct matrices \c Q and \c R, which are resized to the
// correct dimensions (if possible and necessary).
//
// Example:

   \code
   blaze::DynamicMatrix<double,blaze::columnMajor> A( 32, 16 );
   // ... Initialization of A

   blaze::DynamicMatrix<double,blaze::columnMajor> Q( 32, 32 );
   blaze::DynamicMatrix<double,blaze::columnMajor> R( 32, 16 );

   qr( A, Q, R );

   assert( A == Q * R );
   \endcode

// \note This function only works for matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with matrices of any other
// element type results in a compile time error!
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT1  // Type of matrix A
        , bool SO1      // Storage order of matrix A
        , typename MT2  // Type of matrix Q
        , bool SO2      // Storage order of matrix Q
        , typename MT3  // Type of matrix R
        , bool SO3 >    // Storage order of matrix R
void qr( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& Q, DenseMatrix<MT3,SO3>& R )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT1::ElementType );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT2::ElementType );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT3::ElementType );

   typedef typename MT2::ElementType  ET2;
   typedef typename MT3::ElementType  ET3;
   typedef typename RemoveAdaptor<MT3>::Type  UMT3;
   typedef typename If< IsRowMajorMatrix<UMT3>, typename UMT3::OppositeType, UMT3 >::Type  Tmp;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( Tmp );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET3, typename Tmp::ElementType );

   const size_t m( (~A).rows() );
   const size_t n( (~A).columns() );
   const size_t mindim( min( m, n ) );

   if( ( !IsResizable<MT2>::value && ( (~Q).rows() != m || (~Q).columns() != m ) ) ||
       ( !IsResizable<MT3>::value && ( (~R).rows() != m || (~R).columns() != n ) ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Dimensions of fixed size matrix do not match" );
   }

   Tmp tmp( ~A );
   UniqueArray<ET3> tau( new ET3[mindim] );

   geqrf( tmp, tau.get() );

   typename DerestrictTrait<MT3>::Type r( derestrict( ~R ) );
   resize( ~R, m, n );
   reset( r );

   for( size_t i=0UL; i<mindim; ++i ) {
      r(i,i) = tmp(i,i);
      tmp(i,i) = ET3(1);
      for( size_t j=i+1UL; j<n; ++j ) {
         r(i,j) = tmp(i,j);
         reset( tmp(i,j) );
      }
   }

   MT2 I, Qtmp;
   resize( I, m, m );
   reset( I );

   for( size_t i=0UL; i<m; ++i ) {
      I(i,i) = ET2(1);
   }

   (~Q) = I;

   for( size_t i=0UL; i<mindim; ++i ) {
      if( !isDefault( tau[i] ) ) {
         Qtmp = Q;
         DenseColumn<Tmp> col = column( tmp, i );
         (~Q) = Qtmp * ( I - tau[i] * col * ctrans( col ) );
      }
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
