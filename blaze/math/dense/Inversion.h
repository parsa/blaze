//=================================================================================================
/*!
//  \file blaze/math/dense/Inversion.h
//  \brief Header file for the dense matrix inversion functionality
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

#ifndef _BLAZE_MATH_DENSE_INVERSION_H_
#define _BLAZE_MATH_DENSE_INVERSION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/StrictlyTriangular.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/Functions.h>
#include <blaze/math/lapack/Inversion.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsUniTriangular.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Exception.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsFloat.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/UniqueArray.h>


namespace blaze {

//=================================================================================================
//
//  INVERSION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Inversion functions */
//@{
template< typename MT, bool SO >
inline void invert( DenseMatrix<MT,SO>& dm );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense matrix with single precision elements.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the given dense single precision matrix by means of LAPACK kernels. The
// matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \c dm may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename EnableIf< IsFloat<typename MT::ElementType> >::Type
   invert_backend( DenseMatrix<MT,SO>& dm )
{
   const size_t N( min( (~dm).rows(), (~dm).columns() ) );
   UniqueArray<int> ipiv( new int[N] );

   if( IsUniTriangular<MT>::value ) {
      for( size_t i=0UL; i<N; ++i )
         ipiv[i] = static_cast<int>( i ) + 1;
   }
   else {
      sgetrf( derestrict( ~dm ), ipiv.get() );
   }

   sgetri( derestrict( ~dm ), ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense matrix with double precision elements.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the given dense double precision matrix by means of LAPACK kernels.  The
// matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \c dm may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename EnableIf< IsDouble<typename MT::ElementType> >::Type
   invert_backend( DenseMatrix<MT,SO>& dm )
{
   const size_t N( min( (~dm).rows(), (~dm).columns() ) );
   UniqueArray<int> ipiv( new int[N] );

   if( IsUniTriangular<MT>::value ) {
      for( size_t i=0UL; i<N; ++i )
         ipiv[i] = static_cast<int>( i ) + 1;
   }
   else {
      dgetrf( derestrict( ~dm ), ipiv.get() );
   }

   dgetri( derestrict( ~dm ), ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense matrix with single precision complex elements.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the given dense single precision complex matrix by means of LAPACK
// kernels.  The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \c dm may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename EnableIf< IsSame< typename MT::ElementType, complex<float> > >::Type
   invert_backend( DenseMatrix<MT,SO>& dm )
{
   const size_t N( min( (~dm).rows(), (~dm).columns() ) );
   UniqueArray<int> ipiv( new int[N] );

   if( IsUniTriangular<MT>::value ) {
      for( size_t i=0UL; i<N; ++i )
         ipiv[i] = static_cast<int>( i ) + 1;
   }
   else {
      cgetrf( derestrict( ~dm ), ipiv.get() );
   }

   cgetri( derestrict( ~dm ), ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense matrix with double precision complex elements.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the given dense double precision complex matrix by means of LAPACK
// kernels. The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \c dm may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename EnableIf< IsSame< typename MT::ElementType, complex<double> > >::Type
   invert_backend( DenseMatrix<MT,SO>& dm )
{
   const size_t N( min( (~dm).rows(), (~dm).columns() ) );
   UniqueArray<int> ipiv( new int[N] );

   if( IsUniTriangular<MT>::value ) {
      for( size_t i=0UL; i<N; ++i )
         ipiv[i] = static_cast<int>( i ) + 1;
   }
   else {
      zgetrf( derestrict( ~dm ), ipiv.get() );
   }

   zgetri( derestrict( ~dm ), ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the given dense matrix by means of LAPACK kernels. The matrix inversion
// fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or a \a std::invalid_argument exception is thrown.
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \c dm may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invert( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );

   if( !IsSquare<MT>::value && !isSquare( ~dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   invert_backend( ~dm );

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
};
//*************************************************************************************************

} // namespace blaze

#endif
