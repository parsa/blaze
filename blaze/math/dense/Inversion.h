//=================================================================================================
/*!
//  \file blaze/math/dense/Inversion.h
//  \brief Header file for the dense matrix in-place inversion kernels
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

#include <blaze/math/constraints/BlasCompatible.h>
#include <blaze/math/constraints/StrictlyTriangular.h>
#include <blaze/math/DecompositionFlag.h>
#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/Functions.h>
#include <blaze/math/lapack/Cholesky.h>
#include <blaze/math/lapack/Inversion.h>
#include <blaze/math/lapack/PLU.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniTriangular.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Exception.h>
#include <blaze/util/Types.h>
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
// the PLU decomposition or the Cholesky decomposition.
*/
template< DecompositionFlag DF >  // Decomposition algorithm
struct InvertHelper;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DENSE MATRIX INVERSION BASED ON THE PLU DECOMPOSITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the InvertHelper class template for PLU decompositions.
// \ingroup dense_matrix
//
// This specialization of the InvertHelper class template implements the mechanics of the dense
// matrix inversion by means of the PLU decomposition.
*/
template<>
struct InvertHelper<byPLU>
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
// This function inverts the given dense matrix by means of a PLU decomposition. The matrix
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
inline void InvertHelper<byPLU>::invert( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

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
/*!\brief Specialization of the InvertHelper class template for Cholesky decompositions.
// \ingroup dense_matrix
//
// This specialization of the InvertHelper class template implements the mechanics of the dense
// matrix inversion by means of the Cholesky decomposition.
*/
template<>
struct InvertHelper<byCholesky>
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
inline void InvertHelper<byCholesky>::invert( DenseMatrix<MT,SO>& dm )
{
   using blaze::invert;

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( typename MT::ElementType );

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
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
//
// This function inverts the given general dense \f$ 2 \times 2 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invert2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 2UL, "Invalid number of columns detected" );

   typedef typename MT::ElementType  ET;

   typename DerestrictTrait<MT>::Type A( derestrict( ~dm ) );

   const ET det( A(0,0)*A(1,1) - A(0,1)*A(1,0) );

   if( isDefault( det ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );
   const ET a11( A(0,0) * idet );

   A(0,0) =  A(1,1) * idet;
   A(1,0) = -A(1,0) * idet;
   A(0,1) = -A(0,1) * idet;
   A(1,1) =  a11;

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
//
// This function inverts the given general dense \f$ 3 \times 3 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invert3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 3UL, "Invalid number of columns detected" );

   typedef typename MT::ElementType  ET;

   const StaticMatrix<ET,3UL,3UL,SO> A( ~dm );
   typename DerestrictTrait<MT>::Type B( derestrict( ~dm ) );

   B(0,0) = A(1,1)*A(2,2) - A(1,2)*A(2,1);
   B(1,0) = A(1,2)*A(2,0) - A(1,0)*A(2,2);
   B(2,0) = A(1,0)*A(2,1) - A(1,1)*A(2,0);

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) );

   if( isDefault( det ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }

   B(0,1) = A(0,2)*A(2,1) - A(0,1)*A(2,2);
   B(1,1) = A(0,0)*A(2,2) - A(0,2)*A(2,0);
   B(2,1) = A(0,1)*A(2,0) - A(0,0)*A(2,1);
   B(0,2) = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   B(1,2) = A(0,2)*A(1,0) - A(0,0)*A(1,2);
   B(2,2) = A(0,0)*A(1,1) - A(0,1)*A(1,0);

   B /= det;

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
//
// This function inverts the given general dense \f$ 4 \times 4 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invert4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 4UL, "Invalid number of columns detected" );

   typedef typename MT::ElementType  ET;

   const StaticMatrix<ET,4UL,4UL,SO> A( ~dm );
   typename DerestrictTrait<MT>::Type B( derestrict( ~dm ) );

   ET tmp1( A(2,2)*A(3,3) - A(2,3)*A(3,2) );
   ET tmp2( A(2,1)*A(3,3) - A(2,3)*A(3,1) );
   ET tmp3( A(2,1)*A(3,2) - A(2,2)*A(3,1) );
   ET tmp4( A(2,0)*A(3,3) - A(2,3)*A(3,0) );
   ET tmp5( A(2,0)*A(3,2) - A(2,2)*A(3,0) );
   ET tmp6( A(2,0)*A(3,1) - A(2,1)*A(3,0) );

   B(0,0) = A(1,1)*tmp1 - A(1,2)*tmp2 + A(1,3)*tmp3;
   B(1,0) = A(1,2)*tmp4 - A(1,0)*tmp1 - A(1,3)*tmp5;
   B(2,0) = A(1,0)*tmp2 - A(1,1)*tmp4 + A(1,3)*tmp6;
   B(3,0) = A(1,1)*tmp5 - A(1,0)*tmp3 - A(1,2)*tmp6;
   B(0,1) = A(0,2)*tmp2 - A(0,1)*tmp1 - A(0,3)*tmp3;
   B(1,1) = A(0,0)*tmp1 - A(0,2)*tmp4 + A(0,3)*tmp5;
   B(2,1) = A(0,1)*tmp4 - A(0,0)*tmp2 - A(0,3)*tmp6;
   B(3,1) = A(0,0)*tmp3 - A(0,1)*tmp5 + A(0,2)*tmp6;

   tmp1 = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp2 = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp3 = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp4 = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp5 = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp6 = A(0,0)*A(1,1) - A(0,1)*A(1,0);

   B(0,2) = A(3,1)*tmp1 - A(3,2)*tmp2 + A(3,3)*tmp3;
   B(1,2) = A(3,2)*tmp4 - A(3,0)*tmp1 - A(3,3)*tmp5;
   B(2,2) = A(3,0)*tmp2 - A(3,1)*tmp4 + A(3,3)*tmp6;
   B(3,2) = A(3,1)*tmp5 - A(3,0)*tmp3 - A(3,2)*tmp6;
   B(0,3) = A(2,2)*tmp2 - A(2,1)*tmp1 - A(2,3)*tmp3;
   B(1,3) = A(2,0)*tmp1 - A(2,2)*tmp4 + A(2,3)*tmp5;
   B(2,3) = A(2,1)*tmp4 - A(2,0)*tmp2 - A(2,3)*tmp6;
   B(3,3) = A(2,0)*tmp3 - A(2,1)*tmp5 + A(2,2)*tmp6;

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) + A(0,3)*B(3,0) );

   if( isDefault( det ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }

   B /= det;

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
//
// This function inverts the given general dense \f$ 5 \times 5 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invert5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 5UL, "Invalid number of columns detected" );

   typedef typename MT::ElementType  ET;

   const StaticMatrix<ET,5UL,5UL,SO> A( ~dm );
   typename DerestrictTrait<MT>::Type B( derestrict( ~dm ) );

   ET tmp1 ( A(3,3)*A(4,4) - A(3,4)*A(4,3) );
   ET tmp2 ( A(3,2)*A(4,4) - A(3,4)*A(4,2) );
   ET tmp3 ( A(3,2)*A(4,3) - A(3,3)*A(4,2) );
   ET tmp4 ( A(3,1)*A(4,4) - A(3,4)*A(4,1) );
   ET tmp5 ( A(3,1)*A(4,3) - A(3,3)*A(4,1) );
   ET tmp6 ( A(3,1)*A(4,2) - A(3,2)*A(4,1) );
   ET tmp7 ( A(3,0)*A(4,4) - A(3,4)*A(4,0) );
   ET tmp8 ( A(3,0)*A(4,3) - A(3,3)*A(4,0) );
   ET tmp9 ( A(3,0)*A(4,2) - A(3,2)*A(4,0) );
   ET tmp10( A(3,0)*A(4,1) - A(3,1)*A(4,0) );

   ET tmp11( A(2,2)*tmp1 - A(2,3)*tmp2 + A(2,4)*tmp3  );
   ET tmp12( A(2,1)*tmp1 - A(2,3)*tmp4 + A(2,4)*tmp5  );
   ET tmp13( A(2,1)*tmp2 - A(2,2)*tmp4 + A(2,4)*tmp6  );
   ET tmp14( A(2,1)*tmp3 - A(2,2)*tmp5 + A(2,3)*tmp6  );
   ET tmp15( A(1,2)*tmp1 - A(1,3)*tmp2 + A(1,4)*tmp3  );
   ET tmp16( A(1,1)*tmp1 - A(1,3)*tmp4 + A(1,4)*tmp5  );
   ET tmp17( A(1,1)*tmp2 - A(1,2)*tmp4 + A(1,4)*tmp6  );
   ET tmp18( A(1,1)*tmp3 - A(1,2)*tmp5 + A(1,3)*tmp6  );
   ET tmp19( A(2,0)*tmp1 - A(2,3)*tmp7 + A(2,4)*tmp8  );
   ET tmp20( A(2,0)*tmp2 - A(2,2)*tmp7 + A(2,4)*tmp9  );
   ET tmp21( A(2,0)*tmp3 - A(2,2)*tmp8 + A(2,3)*tmp9  );
   ET tmp22( A(1,0)*tmp1 - A(1,3)*tmp7 + A(1,4)*tmp8  );
   ET tmp23( A(1,0)*tmp2 - A(1,2)*tmp7 + A(1,4)*tmp9  );
   ET tmp24( A(1,0)*tmp3 - A(1,2)*tmp8 + A(1,3)*tmp9  );
   ET tmp25( A(2,0)*tmp4 - A(2,1)*tmp7 + A(2,4)*tmp10 );
   ET tmp26( A(2,0)*tmp5 - A(2,1)*tmp8 + A(2,3)*tmp10 );
   ET tmp27( A(1,0)*tmp4 - A(1,1)*tmp7 + A(1,4)*tmp10 );
   ET tmp28( A(1,0)*tmp5 - A(1,1)*tmp8 + A(1,3)*tmp10 );
   ET tmp29( A(2,0)*tmp6 - A(2,1)*tmp9 + A(2,2)*tmp10 );

   B(0,0) =   A(1,1)*tmp11 - A(1,2)*tmp12 + A(1,3)*tmp13 - A(1,4)*tmp14;
   B(1,0) = - A(1,0)*tmp11 + A(1,2)*tmp19 - A(1,3)*tmp20 + A(1,4)*tmp21;
   B(2,0) =   A(1,0)*tmp12 - A(1,1)*tmp19 + A(1,3)*tmp25 - A(1,4)*tmp26;
   B(3,0) = - A(1,0)*tmp13 + A(1,1)*tmp20 - A(1,2)*tmp25 + A(1,4)*tmp29;
   B(4,0) =   A(1,0)*tmp14 - A(1,1)*tmp21 + A(1,2)*tmp26 - A(1,3)*tmp29;
   B(0,1) = - A(0,1)*tmp11 + A(0,2)*tmp12 - A(0,3)*tmp13 + A(0,4)*tmp14;
   B(1,1) =   A(0,0)*tmp11 - A(0,2)*tmp19 + A(0,3)*tmp20 - A(0,4)*tmp21;
   B(2,1) = - A(0,0)*tmp12 + A(0,1)*tmp19 - A(0,3)*tmp25 + A(0,4)*tmp26;
   B(3,1) =   A(0,0)*tmp13 - A(0,1)*tmp20 + A(0,2)*tmp25 - A(0,4)*tmp29;
   B(4,1) = - A(0,0)*tmp14 + A(0,1)*tmp21 - A(0,2)*tmp26 + A(0,3)*tmp29;
   B(0,2) =   A(0,1)*tmp15 - A(0,2)*tmp16 + A(0,3)*tmp17 - A(0,4)*tmp18;
   B(1,2) = - A(0,0)*tmp15 + A(0,2)*tmp22 - A(0,3)*tmp23 + A(0,4)*tmp24;
   B(2,2) =   A(0,0)*tmp16 - A(0,1)*tmp22 + A(0,3)*tmp27 - A(0,4)*tmp28;

   tmp1  = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp2  = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp3  = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp4  = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp5  = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp6  = A(0,0)*A(1,1) - A(0,1)*A(1,0);
   tmp7  = A(0,2)*A(1,4) - A(0,4)*A(1,2);
   tmp8  = A(0,1)*A(1,4) - A(0,4)*A(1,1);
   tmp9  = A(0,0)*A(1,4) - A(0,4)*A(1,0);
   tmp10 = A(0,3)*A(1,4) - A(0,4)*A(1,3);

   tmp11 = A(2,2)*tmp10 - A(2,3)*tmp7 + A(2,4)*tmp1;
   tmp12 = A(2,1)*tmp10 - A(2,3)*tmp8 + A(2,4)*tmp2;
   tmp13 = A(2,1)*tmp7  - A(2,2)*tmp8 + A(2,4)*tmp3;
   tmp14 = A(2,1)*tmp1  - A(2,2)*tmp2 + A(2,3)*tmp3;
   tmp15 = A(2,0)*tmp10 - A(2,3)*tmp9 + A(2,4)*tmp4;
   tmp16 = A(2,0)*tmp7  - A(2,2)*tmp9 + A(2,4)*tmp5;
   tmp17 = A(2,0)*tmp1  - A(2,2)*tmp4 + A(2,3)*tmp5;
   tmp18 = A(2,0)*tmp8  - A(2,1)*tmp9 + A(2,4)*tmp6;
   tmp19 = A(2,0)*tmp2  - A(2,1)*tmp4 + A(2,3)*tmp6;
   tmp20 = A(3,1)*tmp7  - A(3,2)*tmp8 + A(3,4)*tmp3;
   tmp21 = A(3,0)*tmp7  - A(3,2)*tmp9 + A(3,4)*tmp5;
   tmp22 = A(3,0)*tmp8  - A(3,1)*tmp9 + A(3,4)*tmp6;
   tmp23 = A(3,0)*tmp3  - A(3,1)*tmp5 + A(3,2)*tmp6;
   tmp24 = A(2,0)*tmp3  - A(2,1)*tmp5 + A(2,2)*tmp6;
   tmp25 = A(3,1)*tmp1  - A(3,2)*tmp2 + A(3,3)*tmp3;
   tmp26 = A(3,0)*tmp1  - A(3,2)*tmp4 + A(3,3)*tmp5;
   tmp27 = A(3,0)*tmp2  - A(3,1)*tmp4 + A(3,3)*tmp6;

   B(3,2) =   A(4,0)*tmp20 - A(4,1)*tmp21 + A(4,2)*tmp22 - A(4,4)*tmp23;
   B(4,2) = - A(4,0)*tmp25 + A(4,1)*tmp26 - A(4,2)*tmp27 + A(4,3)*tmp23;
   B(0,3) =   A(4,1)*tmp11 - A(4,2)*tmp12 + A(4,3)*tmp13 - A(4,4)*tmp14;
   B(1,3) = - A(4,0)*tmp11 + A(4,2)*tmp15 - A(4,3)*tmp16 + A(4,4)*tmp17;
   B(2,3) =   A(4,0)*tmp12 - A(4,1)*tmp15 + A(4,3)*tmp18 - A(4,4)*tmp19;
   B(3,3) = - A(4,0)*tmp13 + A(4,1)*tmp16 - A(4,2)*tmp18 + A(4,4)*tmp24;
   B(4,3) =   A(4,0)*tmp14 - A(4,1)*tmp17 + A(4,2)*tmp19 - A(4,3)*tmp24;
   B(0,4) = - A(3,1)*tmp11 + A(3,2)*tmp12 - A(3,3)*tmp13 + A(3,4)*tmp14;
   B(1,4) =   A(3,0)*tmp11 - A(3,2)*tmp15 + A(3,3)*tmp16 - A(3,4)*tmp17;
   B(2,4) = - A(3,0)*tmp12 + A(3,1)*tmp15 - A(3,3)*tmp18 + A(3,4)*tmp19;
   B(3,4) =   A(3,0)*tmp13 - A(3,1)*tmp16 + A(3,2)*tmp18 - A(3,4)*tmp24;
   B(4,4) = - A(3,0)*tmp14 + A(3,1)*tmp17 - A(3,2)*tmp19 + A(3,3)*tmp24;

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) + A(0,3)*B(3,0) + A(0,4)*B(4,0) );

   if( isDefault( det ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }

   B /= det;

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
//
// This function inverts the given general dense \f$ 6 \times 6 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invert6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 6UL, "Invalid number of columns detected" );

   typedef typename MT::ElementType  ET;

   const StaticMatrix<ET,6UL,6UL,SO> A( ~dm );
   typename DerestrictTrait<MT>::Type B( derestrict( ~dm ) );

   ET tmp1 ( A(4,4)*A(5,5) - A(4,5)*A(5,4) );
   ET tmp2 ( A(4,3)*A(5,5) - A(4,5)*A(5,3) );
   ET tmp3 ( A(4,3)*A(5,4) - A(4,4)*A(5,3) );
   ET tmp4 ( A(4,2)*A(5,5) - A(4,5)*A(5,2) );
   ET tmp5 ( A(4,2)*A(5,4) - A(4,4)*A(5,2) );
   ET tmp6 ( A(4,2)*A(5,3) - A(4,3)*A(5,2) );
   ET tmp7 ( A(4,1)*A(5,5) - A(4,5)*A(5,1) );
   ET tmp8 ( A(4,1)*A(5,4) - A(4,4)*A(5,1) );
   ET tmp9 ( A(4,1)*A(5,3) - A(4,3)*A(5,1) );
   ET tmp10( A(4,1)*A(5,2) - A(4,2)*A(5,1) );
   ET tmp11( A(4,0)*A(5,5) - A(4,5)*A(5,0) );
   ET tmp12( A(4,0)*A(5,4) - A(4,4)*A(5,0) );
   ET tmp13( A(4,0)*A(5,3) - A(4,3)*A(5,0) );
   ET tmp14( A(4,0)*A(5,2) - A(4,2)*A(5,0) );
   ET tmp15( A(4,0)*A(5,1) - A(4,1)*A(5,0) );

   ET tmp16( A(3,3)*tmp1  - A(3,4)*tmp2  + A(3,5)*tmp3  );
   ET tmp17( A(3,2)*tmp1  - A(3,4)*tmp4  + A(3,5)*tmp5  );
   ET tmp18( A(3,2)*tmp2  - A(3,3)*tmp4  + A(3,5)*tmp6  );
   ET tmp19( A(3,2)*tmp3  - A(3,3)*tmp5  + A(3,4)*tmp6  );
   ET tmp20( A(3,1)*tmp1  - A(3,4)*tmp7  + A(3,5)*tmp8  );
   ET tmp21( A(3,1)*tmp2  - A(3,3)*tmp7  + A(3,5)*tmp9  );
   ET tmp22( A(3,1)*tmp3  - A(3,3)*tmp8  + A(3,4)*tmp9  );
   ET tmp23( A(3,1)*tmp4  - A(3,2)*tmp7  + A(3,5)*tmp10 );
   ET tmp24( A(3,1)*tmp5  - A(3,2)*tmp8  + A(3,4)*tmp10 );
   ET tmp25( A(3,1)*tmp6  - A(3,2)*tmp9  + A(3,3)*tmp10 );
   ET tmp26( A(3,0)*tmp1  - A(3,4)*tmp11 + A(3,5)*tmp12 );
   ET tmp27( A(3,0)*tmp2  - A(3,3)*tmp11 + A(3,5)*tmp13 );
   ET tmp28( A(3,0)*tmp3  - A(3,3)*tmp12 + A(3,4)*tmp13 );
   ET tmp29( A(3,0)*tmp4  - A(3,2)*tmp11 + A(3,5)*tmp14 );
   ET tmp30( A(3,0)*tmp5  - A(3,2)*tmp12 + A(3,4)*tmp14 );
   ET tmp31( A(3,0)*tmp6  - A(3,2)*tmp13 + A(3,3)*tmp14 );
   ET tmp32( A(3,0)*tmp7  - A(3,1)*tmp11 + A(3,5)*tmp15 );
   ET tmp33( A(3,0)*tmp8  - A(3,1)*tmp12 + A(3,4)*tmp15 );
   ET tmp34( A(3,0)*tmp9  - A(3,1)*tmp13 + A(3,3)*tmp15 );
   ET tmp35( A(3,0)*tmp10 - A(3,1)*tmp14 + A(3,2)*tmp15 );

   ET tmp36( A(2,2)*tmp16 - A(2,3)*tmp17 + A(2,4)*tmp18 - A(2,5)*tmp19 );
   ET tmp37( A(2,1)*tmp16 - A(2,3)*tmp20 + A(2,4)*tmp21 - A(2,5)*tmp22 );
   ET tmp38( A(2,1)*tmp17 - A(2,2)*tmp20 + A(2,4)*tmp23 - A(2,5)*tmp24 );
   ET tmp39( A(2,1)*tmp18 - A(2,2)*tmp21 + A(2,3)*tmp23 - A(2,5)*tmp25 );
   ET tmp40( A(2,1)*tmp19 - A(2,2)*tmp22 + A(2,3)*tmp24 - A(2,4)*tmp25 );
   ET tmp41( A(1,2)*tmp16 - A(1,3)*tmp17 + A(1,4)*tmp18 - A(1,5)*tmp19 );
   ET tmp42( A(1,1)*tmp16 - A(1,3)*tmp20 + A(1,4)*tmp21 - A(1,5)*tmp22 );
   ET tmp43( A(1,1)*tmp17 - A(1,2)*tmp20 + A(1,4)*tmp23 - A(1,5)*tmp24 );
   ET tmp44( A(1,1)*tmp18 - A(1,2)*tmp21 + A(1,3)*tmp23 - A(1,5)*tmp25 );
   ET tmp45( A(1,1)*tmp19 - A(1,2)*tmp22 + A(1,3)*tmp24 - A(1,4)*tmp25 );
   ET tmp46( A(2,0)*tmp16 - A(2,3)*tmp26 + A(2,4)*tmp27 - A(2,5)*tmp28 );
   ET tmp47( A(2,0)*tmp17 - A(2,2)*tmp26 + A(2,4)*tmp29 - A(2,5)*tmp30 );
   ET tmp48( A(2,0)*tmp18 - A(2,2)*tmp27 + A(2,3)*tmp29 - A(2,5)*tmp31 );
   ET tmp49( A(2,0)*tmp19 - A(2,2)*tmp28 + A(2,3)*tmp30 - A(2,4)*tmp31 );
   ET tmp50( A(1,0)*tmp16 - A(1,3)*tmp26 + A(1,4)*tmp27 - A(1,5)*tmp28 );
   ET tmp51( A(1,0)*tmp17 - A(1,2)*tmp26 + A(1,4)*tmp29 - A(1,5)*tmp30 );
   ET tmp52( A(1,0)*tmp18 - A(1,2)*tmp27 + A(1,3)*tmp29 - A(1,5)*tmp31 );
   ET tmp53( A(1,0)*tmp19 - A(1,2)*tmp28 + A(1,3)*tmp30 - A(1,4)*tmp31 );
   ET tmp54( A(2,0)*tmp20 - A(2,1)*tmp26 + A(2,4)*tmp32 - A(2,5)*tmp33 );
   ET tmp55( A(2,0)*tmp21 - A(2,1)*tmp27 + A(2,3)*tmp32 - A(2,5)*tmp34 );
   ET tmp56( A(2,0)*tmp22 - A(2,1)*tmp28 + A(2,3)*tmp33 - A(2,4)*tmp34 );
   ET tmp57( A(1,0)*tmp20 - A(1,1)*tmp26 + A(1,4)*tmp32 - A(1,5)*tmp33 );
   ET tmp58( A(1,0)*tmp21 - A(1,1)*tmp27 + A(1,3)*tmp32 - A(1,5)*tmp34 );
   ET tmp59( A(1,0)*tmp22 - A(1,1)*tmp28 + A(1,3)*tmp33 - A(1,4)*tmp34 );
   ET tmp60( A(2,0)*tmp23 - A(2,1)*tmp29 + A(2,2)*tmp32 - A(2,5)*tmp35 );
   ET tmp61( A(2,0)*tmp24 - A(2,1)*tmp30 + A(2,2)*tmp33 - A(2,4)*tmp35 );
   ET tmp62( A(1,0)*tmp23 - A(1,1)*tmp29 + A(1,2)*tmp32 - A(1,5)*tmp35 );
   ET tmp63( A(1,0)*tmp24 - A(1,1)*tmp30 + A(1,2)*tmp33 - A(1,4)*tmp35 );
   ET tmp64( A(2,0)*tmp25 - A(2,1)*tmp31 + A(2,2)*tmp34 - A(2,3)*tmp35 );
   ET tmp65( A(1,0)*tmp25 - A(1,1)*tmp31 + A(1,2)*tmp34 - A(1,3)*tmp35 );

   B(0,0) =   A(1,1)*tmp36 - A(1,2)*tmp37 + A(1,3)*tmp38 - A(1,4)*tmp39 + A(1,5)*tmp40;
   B(1,0) = - A(1,0)*tmp36 + A(1,2)*tmp46 - A(1,3)*tmp47 + A(1,4)*tmp48 - A(1,5)*tmp49;
   B(2,0) =   A(1,0)*tmp37 - A(1,1)*tmp46 + A(1,3)*tmp54 - A(1,4)*tmp55 + A(1,5)*tmp56;
   B(3,0) = - A(1,0)*tmp38 + A(1,1)*tmp47 - A(1,2)*tmp54 + A(1,4)*tmp60 - A(1,5)*tmp61;
   B(4,0) =   A(1,0)*tmp39 - A(1,1)*tmp48 + A(1,2)*tmp55 - A(1,3)*tmp60 + A(1,5)*tmp64;
   B(5,0) = - A(1,0)*tmp40 + A(1,1)*tmp49 - A(1,2)*tmp56 + A(1,3)*tmp61 - A(1,4)*tmp64;
   B(0,1) = - A(0,1)*tmp36 + A(0,2)*tmp37 - A(0,3)*tmp38 + A(0,4)*tmp39 - A(0,5)*tmp40;
   B(1,1) =   A(0,0)*tmp36 - A(0,2)*tmp46 + A(0,3)*tmp47 - A(0,4)*tmp48 + A(0,5)*tmp49;
   B(2,1) = - A(0,0)*tmp37 + A(0,1)*tmp46 - A(0,3)*tmp54 + A(0,4)*tmp55 - A(0,5)*tmp56;
   B(3,1) =   A(0,0)*tmp38 - A(0,1)*tmp47 + A(0,2)*tmp54 - A(0,4)*tmp60 + A(0,5)*tmp61;
   B(4,1) = - A(0,0)*tmp39 + A(0,1)*tmp48 - A(0,2)*tmp55 + A(0,3)*tmp60 - A(0,5)*tmp64;
   B(5,1) =   A(0,0)*tmp40 - A(0,1)*tmp49 + A(0,2)*tmp56 - A(0,3)*tmp61 + A(0,4)*tmp64;
   B(0,2) =   A(0,1)*tmp41 - A(0,2)*tmp42 + A(0,3)*tmp43 - A(0,4)*tmp44 + A(0,5)*tmp45;
   B(1,2) = - A(0,0)*tmp41 + A(0,2)*tmp50 - A(0,3)*tmp51 + A(0,4)*tmp52 - A(0,5)*tmp53;
   B(2,2) =   A(0,0)*tmp42 - A(0,1)*tmp50 + A(0,3)*tmp57 - A(0,4)*tmp58 + A(0,5)*tmp59;
   B(3,2) = - A(0,0)*tmp43 + A(0,1)*tmp51 - A(0,2)*tmp57 + A(0,4)*tmp62 - A(0,5)*tmp63;
   B(4,2) =   A(0,0)*tmp44 - A(0,1)*tmp52 + A(0,2)*tmp58 - A(0,3)*tmp62 + A(0,5)*tmp65;
   B(5,2) = - A(0,0)*tmp45 + A(0,1)*tmp53 - A(0,2)*tmp59 + A(0,3)*tmp63 - A(0,4)*tmp65;

   tmp1  = A(0,3)*A(1,4) - A(0,4)*A(1,3);
   tmp2  = A(0,2)*A(1,4) - A(0,4)*A(1,2);
   tmp3  = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp4  = A(0,1)*A(1,4) - A(0,4)*A(1,1);
   tmp5  = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp6  = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp7  = A(0,0)*A(1,4) - A(0,4)*A(1,0);
   tmp8  = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp9  = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp10 = A(0,0)*A(1,1) - A(0,1)*A(1,0);
   tmp11 = A(0,3)*A(1,5) - A(0,5)*A(1,3);
   tmp12 = A(0,2)*A(1,5) - A(0,5)*A(1,2);
   tmp13 = A(0,1)*A(1,5) - A(0,5)*A(1,1);
   tmp14 = A(0,0)*A(1,5) - A(0,5)*A(1,0);
   tmp15 = A(0,4)*A(1,5) - A(0,5)*A(1,4);

   tmp16 = A(2,3)*tmp15 - A(2,4)*tmp11 + A(2,5)*tmp1;
   tmp17 = A(2,2)*tmp15 - A(2,4)*tmp12 + A(2,5)*tmp2;
   tmp18 = A(2,2)*tmp11 - A(2,3)*tmp12 + A(2,5)*tmp3;
   tmp19 = A(2,2)*tmp1  - A(2,3)*tmp2  + A(2,4)*tmp3;
   tmp20 = A(2,1)*tmp15 - A(2,4)*tmp13 + A(2,5)*tmp4;
   tmp21 = A(2,1)*tmp11 - A(2,3)*tmp13 + A(2,5)*tmp5;
   tmp22 = A(2,1)*tmp1  - A(2,3)*tmp4  + A(2,4)*tmp5;
   tmp23 = A(2,1)*tmp12 - A(2,2)*tmp13 + A(2,5)*tmp6;
   tmp24 = A(2,1)*tmp2  - A(2,2)*tmp4  + A(2,4)*tmp6;
   tmp25 = A(2,1)*tmp3  - A(2,2)*tmp5  + A(2,3)*tmp6;
   tmp26 = A(2,0)*tmp15 - A(2,4)*tmp14 + A(2,5)*tmp7;
   tmp27 = A(2,0)*tmp11 - A(2,3)*tmp14 + A(2,5)*tmp8;
   tmp28 = A(2,0)*tmp1  - A(2,3)*tmp7  + A(2,4)*tmp8;
   tmp29 = A(2,0)*tmp12 - A(2,2)*tmp14 + A(2,5)*tmp9;
   tmp30 = A(2,0)*tmp2  - A(2,2)*tmp7  + A(2,4)*tmp9;
   tmp31 = A(2,0)*tmp3  - A(2,2)*tmp8  + A(2,3)*tmp9;
   tmp32 = A(2,0)*tmp13 - A(2,1)*tmp14 + A(2,5)*tmp10;
   tmp33 = A(2,0)*tmp4  - A(2,1)*tmp7  + A(2,4)*tmp10;
   tmp34 = A(2,0)*tmp5  - A(2,1)*tmp8  + A(2,3)*tmp10;
   tmp35 = A(2,0)*tmp6  - A(2,1)*tmp9  + A(2,2)*tmp10;

   tmp36 = A(4,2)*tmp16 - A(4,3)*tmp17 + A(4,4)*tmp18 - A(4,5)*tmp19;
   tmp37 = A(4,1)*tmp16 - A(4,3)*tmp20 + A(4,4)*tmp21 - A(4,5)*tmp22;
   tmp38 = A(4,1)*tmp17 - A(4,2)*tmp20 + A(4,4)*tmp23 - A(4,5)*tmp24;
   tmp39 = A(4,1)*tmp18 - A(4,2)*tmp21 + A(4,3)*tmp23 - A(4,5)*tmp25;
   tmp40 = A(4,1)*tmp19 - A(4,2)*tmp22 + A(4,3)*tmp24 - A(4,4)*tmp25;
   tmp41 = A(3,2)*tmp16 - A(3,3)*tmp17 + A(3,4)*tmp18 - A(3,5)*tmp19;
   tmp42 = A(3,1)*tmp16 - A(3,3)*tmp20 + A(3,4)*tmp21 - A(3,5)*tmp22;
   tmp43 = A(3,1)*tmp17 - A(3,2)*tmp20 + A(3,4)*tmp23 - A(3,5)*tmp24;
   tmp44 = A(3,1)*tmp18 - A(3,2)*tmp21 + A(3,3)*tmp23 - A(3,5)*tmp25;
   tmp45 = A(3,1)*tmp19 - A(3,2)*tmp22 + A(3,3)*tmp24 - A(3,4)*tmp25;
   tmp46 = A(4,0)*tmp16 - A(4,3)*tmp26 + A(4,4)*tmp27 - A(4,5)*tmp28;
   tmp47 = A(4,0)*tmp17 - A(4,2)*tmp26 + A(4,4)*tmp29 - A(4,5)*tmp30;
   tmp48 = A(4,0)*tmp18 - A(4,2)*tmp27 + A(4,3)*tmp29 - A(4,5)*tmp31;
   tmp49 = A(4,0)*tmp19 - A(4,2)*tmp28 + A(4,3)*tmp30 - A(4,4)*tmp31;
   tmp50 = A(3,0)*tmp16 - A(3,3)*tmp26 + A(3,4)*tmp27 - A(3,5)*tmp28;
   tmp51 = A(3,0)*tmp17 - A(3,2)*tmp26 + A(3,4)*tmp29 - A(3,5)*tmp30;
   tmp52 = A(3,0)*tmp18 - A(3,2)*tmp27 + A(3,3)*tmp29 - A(3,5)*tmp31;
   tmp53 = A(3,0)*tmp19 - A(3,2)*tmp28 + A(3,3)*tmp30 - A(3,4)*tmp31;
   tmp54 = A(4,0)*tmp20 - A(4,1)*tmp26 + A(4,4)*tmp32 - A(4,5)*tmp33;
   tmp55 = A(4,0)*tmp21 - A(4,1)*tmp27 + A(4,3)*tmp32 - A(4,5)*tmp34;
   tmp56 = A(4,0)*tmp22 - A(4,1)*tmp28 + A(4,3)*tmp33 - A(4,4)*tmp34;
   tmp57 = A(3,0)*tmp20 - A(3,1)*tmp26 + A(3,4)*tmp32 - A(3,5)*tmp33;
   tmp58 = A(3,0)*tmp21 - A(3,1)*tmp27 + A(3,3)*tmp32 - A(3,5)*tmp34;
   tmp59 = A(3,0)*tmp22 - A(3,1)*tmp28 + A(3,3)*tmp33 - A(3,4)*tmp34;
   tmp60 = A(4,0)*tmp23 - A(4,1)*tmp29 + A(4,2)*tmp32 - A(4,5)*tmp35;
   tmp61 = A(4,0)*tmp24 - A(4,1)*tmp30 + A(4,2)*tmp33 - A(4,4)*tmp35;
   tmp62 = A(3,0)*tmp23 - A(3,1)*tmp29 + A(3,2)*tmp32 - A(3,5)*tmp35;
   tmp63 = A(3,0)*tmp24 - A(3,1)*tmp30 + A(3,2)*tmp33 - A(3,4)*tmp35;
   tmp64 = A(4,0)*tmp25 - A(4,1)*tmp31 + A(4,2)*tmp34 - A(4,3)*tmp35;
   tmp65 = A(3,0)*tmp25 - A(3,1)*tmp31 + A(3,2)*tmp34 - A(3,3)*tmp35;

   B(0,3) =   A(5,1)*tmp36 - A(5,2)*tmp37 + A(5,3)*tmp38 - A(5,4)*tmp39 + A(5,5)*tmp40;
   B(1,3) = - A(5,0)*tmp36 + A(5,2)*tmp46 - A(5,3)*tmp47 + A(5,4)*tmp48 - A(5,5)*tmp49;
   B(2,3) =   A(5,0)*tmp37 - A(5,1)*tmp46 + A(5,3)*tmp54 - A(5,4)*tmp55 + A(5,5)*tmp56;
   B(3,3) = - A(5,0)*tmp38 + A(5,1)*tmp47 - A(5,2)*tmp54 + A(5,4)*tmp60 - A(5,5)*tmp61;
   B(4,3) =   A(5,0)*tmp39 - A(5,1)*tmp48 + A(5,2)*tmp55 - A(5,3)*tmp60 + A(5,5)*tmp64;
   B(5,3) = - A(5,0)*tmp40 + A(5,1)*tmp49 - A(5,2)*tmp56 + A(5,3)*tmp61 - A(5,4)*tmp64;
   B(0,4) = - A(5,1)*tmp41 + A(5,2)*tmp42 - A(5,3)*tmp43 + A(5,4)*tmp44 - A(5,5)*tmp45;
   B(1,4) =   A(5,0)*tmp41 - A(5,2)*tmp50 + A(5,3)*tmp51 - A(5,4)*tmp52 + A(5,5)*tmp53;
   B(2,4) = - A(5,0)*tmp42 + A(5,1)*tmp50 - A(5,3)*tmp57 + A(5,4)*tmp58 - A(5,5)*tmp59;
   B(3,4) =   A(5,0)*tmp43 - A(5,1)*tmp51 + A(5,2)*tmp57 - A(5,4)*tmp62 + A(5,5)*tmp63;
   B(4,4) = - A(5,0)*tmp44 + A(5,1)*tmp52 - A(5,2)*tmp58 + A(5,3)*tmp62 - A(5,5)*tmp65;
   B(5,4) =   A(5,0)*tmp45 - A(5,1)*tmp53 + A(5,2)*tmp59 - A(5,3)*tmp63 + A(5,4)*tmp65;
   B(0,5) =   A(4,1)*tmp41 - A(4,2)*tmp42 + A(4,3)*tmp43 - A(4,4)*tmp44 + A(4,5)*tmp45;
   B(1,5) = - A(4,0)*tmp41 + A(4,2)*tmp50 - A(4,3)*tmp51 + A(4,4)*tmp52 - A(4,5)*tmp53;
   B(2,5) =   A(4,0)*tmp42 - A(4,1)*tmp50 + A(4,3)*tmp57 - A(4,4)*tmp58 + A(4,5)*tmp59;
   B(3,5) = - A(4,0)*tmp43 + A(4,1)*tmp51 - A(4,2)*tmp57 + A(4,4)*tmp62 - A(4,5)*tmp63;
   B(4,5) =   A(4,0)*tmp44 - A(4,1)*tmp52 + A(4,2)*tmp58 - A(4,3)*tmp62 + A(4,5)*tmp65;
   B(5,5) = - A(4,0)*tmp45 + A(4,1)*tmp53 - A(4,2)*tmp59 + A(4,3)*tmp63 - A(4,4)*tmp65;

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) +
                 A(0,3)*B(3,0) + A(0,4)*B(4,0) + A(0,5)*B(5,0) );

   if( isDefault( det ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Inversion of singular matrix failed" );
   }

   B /= det;

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense square matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
//
// This function inverts the given dense square matrix via the specified matrix decomposition
// algorithm \a DF. In case the given matrix is a positive-definite matrix it is recommended
// to perform the inversion by means of a Cholesky decomposition, for a general matrix a PLU
// decomposition should be used:

   \code
   invertNxN<byPLU>( A );       // Inversion of a general matrix
   invertNxN<byCholesky>( A );  // Inversion of a positive definite matrix
   \endcode

// The matrix inversion fails if the given matrix is singular and not invertible. In this case
// a \a std::invalid_argument exception is thrown.
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown, \c m may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< DecompositionFlag DF  // Decomposition algorithm
        , typename MT           // Type of the dense matrix
        , bool SO >             // Storage order of the dense matrix
inline void invertNxN( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( isSquare( ~dm ), "Non-square matrix detected" );

   InvertHelper<DF>::invert( ~dm );

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
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
   invert<byPLU>( ~dm );
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
// a PLU decomposition should be used:

   \code
   invert<byPLU>( A );       // Inversion of a general square matrix
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

   const size_t N( (~dm).rows() );

   if( N == 1UL )
   {
      invert( (~dm)(0,0) );
   }
   else if( N == 2UL )
   {
      invert2x2( ~dm );
   }
   else if( N == 3UL )
   {
      invert3x3( ~dm );
   }
   else if( N == 4UL )
   {
      invert4x4( ~dm );
   }
   else if( N == 5UL )
   {
      invert5x5( ~dm );
   }
   else if( N == 6UL )
   {
      invert6x6( ~dm );
   }
   else
   {
      invertNxN<DF>( ~dm );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( ~dm ), "Broken invariant detected" );
};
//*************************************************************************************************

} // namespace blaze

#endif
