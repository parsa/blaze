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
#include <blaze/math/DecompositionFlag.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/Functions.h>
#include <blaze/math/lapack/Cholesky.h>
#include <blaze/math/lapack/Inversion.h>
#include <blaze/math/lapack/LU.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsTriangular.h>
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
//  CLASS INVERSION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary class template for the implementation of different inversion algorithms.
// \ingroup dense_matrix
//
// This class template represents the base template for the implementation of different dense
// matrix inversion algorithms. In order to implement a specific algorithm this base template
// needs to be specialized for a specific dense matrix decomposition algorithm, as for instance
// the LU decomposition or the Cholesky decomposition.
*/
template< DecompositionFlag DF >  // Decomposition algorithm
struct Inversion;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DENSE MATRIX INVERSION BASED ON THE LU DECOMPOSITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Inversion class template for LU decompositions.
// \ingroup dense_matrix
//
// This specialization of the Inversion class template implements the mechanics of the dense
// matrix inversion by means of the LU decomposition.
*/
template<>
struct Inversion<byLU>
{
   //**Invert functions****************************************************************************
   /*!\name Invert functions */
   //@{
   template< typename MT, bool SO >
   static inline void invert( DenseMatrix<MT,SO>& dm );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the given dense matrix by means of a LU decomposition. The matrix
// inversion fails if ...
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
inline void Inversion<byLU>::invert( DenseMatrix<MT,SO>& dm )
{
   const size_t N( min( (~dm).rows(), (~dm).columns() ) );
   UniqueArray<int> ipiv( new int[N] );

   typename DerestrictTrait<MT>::Type A( derestrict( ~dm ) );

   if( IsUniTriangular<MT>::value ) {
      for( size_t i=0UL; i<N; ++i )
         ipiv[i] = static_cast<int>( i ) + 1;
   }
   else {
      getrf( A, ipiv.get() );
   }

   getri( A, ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DENSE MATRIX INVERSION BASED ON THE CHOLESKY DECOMPOSITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Inversion class template for Cholesky decompositions.
// \ingroup dense_matrix
//
// This specialization of the Inversion class template implements the mechanics of the dense
// matrix inversion by means of the Cholesky decomposition.
*/
template<>
struct Inversion<byCholesky>
{
   //**Invert functions****************************************************************************
   /*!\name Invert functions */
   //@{
   template< typename MT, bool SO >
   static inline void invert( DenseMatrix<MT,SO>& dm );
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the given dense matrix by means of a Cholesky decomposition. The matrix
// inversion fails if ...
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
inline void Inversion<byCholesky>::invert( DenseMatrix<MT,SO>& dm )
{
   using blaze::invert;

   BLAZE_USER_ASSERT( isSymmetric( ~dm ), "Invalid non-symmetric matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( ~dm ), "Non-square matrix detected" );

   if( IsUniTriangular<MT>::value )
      return;

   typename DerestrictTrait<MT>::Type A( derestrict( ~dm ) );

   if( IsTriangular<MT>::value )
   {
      for( size_t i=0UL; i<A.rows(); ++i )
      {
         if( isDefault( A(i,i) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
         }

         invert( A(i,i) );
      }
   }
   else
   {
      const char uplo( ( SO )?( 'L' ):( 'U' ) );

      potrf( A, uplo );
      potri( A, uplo );

      if( SO ) {
         for( size_t i=1UL; i<A.rows(); ++i ) {
            for( size_t j=0UL; j<i; ++j ) {
               A(j,i) = A(i,j);
            }
         }
      }
      else {
         for( size_t j=1UL; j<A.columns(); ++j ) {
            for( size_t i=0UL; i<j; ++i ) {
               A(j,i) = A(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************




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

template< DecompositionFlag DF, typename MT, bool SO >
inline void invert( DenseMatrix<MT,SO>& dm );
//@}
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
// This function inverts the given dense square matrix. The matrix inversion fails if ...
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
   invert<byLU>( ~dm );
};
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
// This function inverts the given dense matrix by means of the specified matrix decomposition
// algorithm \ DF. In case the matrix is a symmetric positive-definite matrix it is recommended
// to perform the inversion by means of a Cholesky decomposition, for a general square matrix
// an LU decomposition should be used:

   \code
   invert<byLU>( A );        // Inversion of a general square matrix
   invert<byCholesky>( A );  // Inversion of a positive definite matrix
   \endcode

// The matrix inversion fails if ...
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
template< DecompositionFlag DF  // Decomposition algorithm
        , typename MT           // Type of the dense matrix
        , bool SO >             // Storage order of the dense matrix
inline void invert( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );

   if( !isSquare( ~dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   Inversion<DF>::invert( ~dm );

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
};
//*************************************************************************************************

} // namespace blaze

#endif
