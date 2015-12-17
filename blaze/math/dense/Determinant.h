//=================================================================================================
/*!
//  \file blaze/math/dense/Determinant.h
//  \brief Header file for the dense matrix determinant functionality
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

#ifndef _BLAZE_MATH_DENSE_DETERMINANT_H_
#define _BLAZE_MATH_DENSE_DETERMINANT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/DecompositionFlag.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/lapack/LU.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsStrictlyTriangular.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniTriangular.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Exception.h>
#include <blaze/util/Types.h>
#include <blaze/util/UniqueArray.h>


namespace blaze {

//=================================================================================================
//
//  DETERMINANT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Determinant functions */
//@{
template< typename MT, bool SO >
inline typename MT::ElementType det( const DenseMatrix<MT,SO>& dm );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The determinant of the given matrix.
//
// This function computes the determinant of the given dense \f$ 2 \times 2 \f$ matrix via the
// rule of Sarrus.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename MT::ElementType det2x2( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 2UL, "Invalid number of columns detected" );

   return (~dm)(0,0)*(~dm)(1,1) - (~dm)(0,1)*(~dm)(1,0);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The determinant of the given matrix.
//
// This function computes the determinant of the given dense \f$ 3 \times 3 \f$ matrix via the
// rule of Sarrus.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename MT::ElementType det3x3( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 3UL, "Invalid number of columns detected" );

   return (~dm)(0,0) * ( (~dm)(1,1)*(~dm)(2,2) - (~dm)(1,2)*(~dm)(2,1) ) +
          (~dm)(0,1) * ( (~dm)(1,2)*(~dm)(2,0) - (~dm)(1,0)*(~dm)(2,2) ) +
          (~dm)(0,2) * ( (~dm)(1,0)*(~dm)(2,1) - (~dm)(1,1)*(~dm)(2,0) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The determinant of the given matrix.
//
// This function computes the determinant of the given dense \f$ 4 \times 4 \f$ matrix via the
// rule of Sarrus.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename MT::ElementType det4x4( const DenseMatrix<MT,SO>& dm )
{
   typedef typename MT::ElementType  ET;

   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 4UL, "Invalid number of columns detected" );

   const ET tmp1( (~dm)(2,2)*(~dm)(3,3) - (~dm)(2,3)*(~dm)(3,2) );
   const ET tmp2( (~dm)(2,1)*(~dm)(3,3) - (~dm)(2,3)*(~dm)(3,1) );
   const ET tmp3( (~dm)(2,1)*(~dm)(3,2) - (~dm)(2,2)*(~dm)(3,1) );
   const ET tmp4( (~dm)(2,0)*(~dm)(3,3) - (~dm)(2,3)*(~dm)(3,0) );
   const ET tmp5( (~dm)(2,0)*(~dm)(3,2) - (~dm)(2,2)*(~dm)(3,0) );
   const ET tmp6( (~dm)(2,0)*(~dm)(3,1) - (~dm)(2,1)*(~dm)(3,0) );

   return (~dm)(0,0) * ( (~dm)(1,1) * tmp1 - (~dm)(1,2) * tmp2 + (~dm)(1,3) * tmp3 ) -
          (~dm)(0,1) * ( (~dm)(1,0) * tmp1 - (~dm)(1,2) * tmp4 + (~dm)(1,3) * tmp5 ) +
          (~dm)(0,2) * ( (~dm)(1,0) * tmp2 - (~dm)(1,1) * tmp4 + (~dm)(1,3) * tmp6 ) -
          (~dm)(0,3) * ( (~dm)(1,0) * tmp3 - (~dm)(1,1) * tmp5 + (~dm)(1,2) * tmp6 );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The determinant of the given matrix.
//
// This function computes the determinant of the given dense \f$ 5 \times 5 \f$ matrix via the
// rule of Sarrus.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename MT::ElementType det5x5( const DenseMatrix<MT,SO>& dm )
{
   typedef typename MT::ElementType  ET;

   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 5UL, "Invalid number of columns detected" );

   const ET tmp1 ( (~dm)(3,3)*(~dm)(4,4) - (~dm)(3,4)*(~dm)(4,3) );
   const ET tmp2 ( (~dm)(3,2)*(~dm)(4,4) - (~dm)(3,4)*(~dm)(4,2) );
   const ET tmp3 ( (~dm)(3,2)*(~dm)(4,3) - (~dm)(3,3)*(~dm)(4,2) );
   const ET tmp4 ( (~dm)(3,1)*(~dm)(4,4) - (~dm)(3,4)*(~dm)(4,1) );
   const ET tmp5 ( (~dm)(3,1)*(~dm)(4,3) - (~dm)(3,3)*(~dm)(4,1) );
   const ET tmp6 ( (~dm)(3,1)*(~dm)(4,2) - (~dm)(3,2)*(~dm)(4,1) );
   const ET tmp7 ( (~dm)(3,0)*(~dm)(4,4) - (~dm)(3,4)*(~dm)(4,0) );
   const ET tmp8 ( (~dm)(3,0)*(~dm)(4,3) - (~dm)(3,3)*(~dm)(4,0) );
   const ET tmp9 ( (~dm)(3,0)*(~dm)(4,2) - (~dm)(3,2)*(~dm)(4,0) );
   const ET tmp10( (~dm)(3,0)*(~dm)(4,1) - (~dm)(3,1)*(~dm)(4,0) );

   const ET tmp11( (~dm)(2,2)*tmp1 - (~dm)(2,3)*tmp2 + (~dm)(2,4)*tmp3 );
   const ET tmp12( (~dm)(2,1)*tmp1 - (~dm)(2,3)*tmp4 + (~dm)(2,4)*tmp5 );
   const ET tmp13( (~dm)(2,1)*tmp2 - (~dm)(2,2)*tmp4 + (~dm)(2,4)*tmp6 );
   const ET tmp14( (~dm)(2,1)*tmp3 - (~dm)(2,2)*tmp5 + (~dm)(2,3)*tmp6 );
   const ET tmp15( (~dm)(2,0)*tmp1 - (~dm)(2,3)*tmp7 + (~dm)(2,4)*tmp8 );
   const ET tmp16( (~dm)(2,0)*tmp2 - (~dm)(2,2)*tmp7 + (~dm)(2,4)*tmp9 );
   const ET tmp17( (~dm)(2,0)*tmp3 - (~dm)(2,2)*tmp8 + (~dm)(2,3)*tmp9 );
   const ET tmp18( (~dm)(2,0)*tmp4 - (~dm)(2,1)*tmp7 + (~dm)(2,4)*tmp10 );
   const ET tmp19( (~dm)(2,0)*tmp5 - (~dm)(2,1)*tmp8 + (~dm)(2,3)*tmp10 );
   const ET tmp20( (~dm)(2,0)*tmp6 - (~dm)(2,1)*tmp9 + (~dm)(2,2)*tmp10 );

   return (~dm)(0,0) * ( (~dm)(1,1)*tmp11 - (~dm)(1,2)*tmp12 + (~dm)(1,3)*tmp13 - (~dm)(1,4)*tmp14 ) -
          (~dm)(0,1) * ( (~dm)(1,0)*tmp11 - (~dm)(1,2)*tmp15 + (~dm)(1,3)*tmp16 - (~dm)(1,4)*tmp17 ) +
          (~dm)(0,2) * ( (~dm)(1,0)*tmp12 - (~dm)(1,1)*tmp15 + (~dm)(1,3)*tmp18 - (~dm)(1,4)*tmp19 ) -
          (~dm)(0,3) * ( (~dm)(1,0)*tmp13 - (~dm)(1,1)*tmp16 + (~dm)(1,2)*tmp18 - (~dm)(1,4)*tmp20 ) +
          (~dm)(0,4) * ( (~dm)(1,0)*tmp14 - (~dm)(1,1)*tmp17 + (~dm)(1,2)*tmp19 - (~dm)(1,3)*tmp20 );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The determinant of the given matrix.
//
// This function computes the determinant of the given dense \f$ 6 \times 6 \f$ matrix via the
// rule of Sarrus.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename MT::ElementType det6x6( const DenseMatrix<MT,SO>& dm )
{
   typedef typename MT::ElementType  ET;

   BLAZE_INTERNAL_ASSERT( (~dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (~dm).columns() == 6UL, "Invalid number of columns detected" );

   const ET tmp1 ( (~dm)(4,4)*(~dm)(5,5) - (~dm)(4,5)*(~dm)(5,4) );
   const ET tmp2 ( (~dm)(4,3)*(~dm)(5,5) - (~dm)(4,5)*(~dm)(5,3) );
   const ET tmp3 ( (~dm)(4,3)*(~dm)(5,4) - (~dm)(4,4)*(~dm)(5,3) );
   const ET tmp4 ( (~dm)(4,2)*(~dm)(5,5) - (~dm)(4,5)*(~dm)(5,2) );
   const ET tmp5 ( (~dm)(4,2)*(~dm)(5,4) - (~dm)(4,4)*(~dm)(5,2) );
   const ET tmp6 ( (~dm)(4,2)*(~dm)(5,3) - (~dm)(4,3)*(~dm)(5,2) );
   const ET tmp7 ( (~dm)(4,1)*(~dm)(5,5) - (~dm)(4,5)*(~dm)(5,1) );
   const ET tmp8 ( (~dm)(4,1)*(~dm)(5,4) - (~dm)(4,4)*(~dm)(5,1) );
   const ET tmp9 ( (~dm)(4,1)*(~dm)(5,3) - (~dm)(4,3)*(~dm)(5,1) );
   const ET tmp10( (~dm)(4,1)*(~dm)(5,2) - (~dm)(4,2)*(~dm)(5,1) );
   const ET tmp11( (~dm)(4,0)*(~dm)(5,5) - (~dm)(4,5)*(~dm)(5,0) );
   const ET tmp12( (~dm)(4,0)*(~dm)(5,4) - (~dm)(4,4)*(~dm)(5,0) );
   const ET tmp13( (~dm)(4,0)*(~dm)(5,3) - (~dm)(4,3)*(~dm)(5,0) );
   const ET tmp14( (~dm)(4,0)*(~dm)(5,2) - (~dm)(4,2)*(~dm)(5,0) );
   const ET tmp15( (~dm)(4,0)*(~dm)(5,1) - (~dm)(4,1)*(~dm)(5,0) );

   const ET tmp16( (~dm)(3,3)*tmp1 - (~dm)(3,4)*tmp2 + (~dm)(3,5)*tmp3 );
   const ET tmp17( (~dm)(3,2)*tmp1 - (~dm)(3,4)*tmp4 + (~dm)(3,5)*tmp5 );
   const ET tmp18( (~dm)(3,2)*tmp2 - (~dm)(3,3)*tmp4 + (~dm)(3,5)*tmp6 );
   const ET tmp19( (~dm)(3,2)*tmp3 - (~dm)(3,3)*tmp5 + (~dm)(3,4)*tmp6 );
   const ET tmp20( (~dm)(3,1)*tmp1 - (~dm)(3,4)*tmp7 + (~dm)(3,5)*tmp8 );
   const ET tmp21( (~dm)(3,1)*tmp2 - (~dm)(3,3)*tmp7 + (~dm)(3,5)*tmp9 );
   const ET tmp22( (~dm)(3,1)*tmp3 - (~dm)(3,3)*tmp8 + (~dm)(3,4)*tmp9 );
   const ET tmp23( (~dm)(3,1)*tmp4 - (~dm)(3,2)*tmp7 + (~dm)(3,5)*tmp10 );
   const ET tmp24( (~dm)(3,1)*tmp5 - (~dm)(3,2)*tmp8 + (~dm)(3,4)*tmp10 );
   const ET tmp25( (~dm)(3,1)*tmp6 - (~dm)(3,2)*tmp9 + (~dm)(3,3)*tmp10 );
   const ET tmp26( (~dm)(3,0)*tmp1 - (~dm)(3,4)*tmp11 + (~dm)(3,5)*tmp12 );
   const ET tmp27( (~dm)(3,0)*tmp2 - (~dm)(3,3)*tmp11 + (~dm)(3,5)*tmp13 );
   const ET tmp28( (~dm)(3,0)*tmp3 - (~dm)(3,3)*tmp12 + (~dm)(3,4)*tmp13 );
   const ET tmp29( (~dm)(3,0)*tmp4 - (~dm)(3,2)*tmp11 + (~dm)(3,5)*tmp14 );
   const ET tmp30( (~dm)(3,0)*tmp5 - (~dm)(3,2)*tmp12 + (~dm)(3,4)*tmp14 );
   const ET tmp31( (~dm)(3,0)*tmp6 - (~dm)(3,2)*tmp13 + (~dm)(3,3)*tmp14 );
   const ET tmp32( (~dm)(3,0)*tmp7 - (~dm)(3,1)*tmp11 + (~dm)(3,5)*tmp15 );
   const ET tmp33( (~dm)(3,0)*tmp8 - (~dm)(3,1)*tmp12 + (~dm)(3,4)*tmp15 );
   const ET tmp34( (~dm)(3,0)*tmp9 - (~dm)(3,1)*tmp13 + (~dm)(3,3)*tmp15 );
   const ET tmp35( (~dm)(3,0)*tmp10 - (~dm)(3,1)*tmp14 + (~dm)(3,2)*tmp15 );

   const ET tmp36( (~dm)(2,2)*tmp16 - (~dm)(2,3)*tmp17 + (~dm)(2,4)*tmp18 - (~dm)(2,5)*tmp19 );
   const ET tmp37( (~dm)(2,1)*tmp16 - (~dm)(2,3)*tmp20 + (~dm)(2,4)*tmp21 - (~dm)(2,5)*tmp22 );
   const ET tmp38( (~dm)(2,1)*tmp17 - (~dm)(2,2)*tmp20 + (~dm)(2,4)*tmp23 - (~dm)(2,5)*tmp24 );
   const ET tmp39( (~dm)(2,1)*tmp18 - (~dm)(2,2)*tmp21 + (~dm)(2,3)*tmp23 - (~dm)(2,5)*tmp25 );
   const ET tmp40( (~dm)(2,1)*tmp19 - (~dm)(2,2)*tmp22 + (~dm)(2,3)*tmp24 - (~dm)(2,4)*tmp25 );
   const ET tmp41( (~dm)(2,0)*tmp16 - (~dm)(2,3)*tmp26 + (~dm)(2,4)*tmp27 - (~dm)(2,5)*tmp28 );
   const ET tmp42( (~dm)(2,0)*tmp17 - (~dm)(2,2)*tmp26 + (~dm)(2,4)*tmp29 - (~dm)(2,5)*tmp30 );
   const ET tmp43( (~dm)(2,0)*tmp18 - (~dm)(2,2)*tmp27 + (~dm)(2,3)*tmp29 - (~dm)(2,5)*tmp31 );
   const ET tmp44( (~dm)(2,0)*tmp19 - (~dm)(2,2)*tmp28 + (~dm)(2,3)*tmp30 - (~dm)(2,4)*tmp31 );
   const ET tmp45( (~dm)(2,0)*tmp20 - (~dm)(2,1)*tmp26 + (~dm)(2,4)*tmp32 - (~dm)(2,5)*tmp33 );
   const ET tmp46( (~dm)(2,0)*tmp21 - (~dm)(2,1)*tmp27 + (~dm)(2,3)*tmp32 - (~dm)(2,5)*tmp34 );
   const ET tmp47( (~dm)(2,0)*tmp22 - (~dm)(2,1)*tmp28 + (~dm)(2,3)*tmp33 - (~dm)(2,4)*tmp34 );
   const ET tmp48( (~dm)(2,0)*tmp23 - (~dm)(2,1)*tmp29 + (~dm)(2,2)*tmp32 - (~dm)(2,5)*tmp35 );
   const ET tmp49( (~dm)(2,0)*tmp24 - (~dm)(2,1)*tmp30 + (~dm)(2,2)*tmp33 - (~dm)(2,4)*tmp35 );
   const ET tmp50( (~dm)(2,0)*tmp25 - (~dm)(2,1)*tmp31 + (~dm)(2,2)*tmp34 - (~dm)(2,3)*tmp35 );

   return (~dm)(0,0) * ( (~dm)(1,1)*tmp36 - (~dm)(1,2)*tmp37 + (~dm)(1,3)*tmp38 - (~dm)(1,4)*tmp39 + (~dm)(1,5)*tmp40 ) -
          (~dm)(0,1) * ( (~dm)(1,0)*tmp36 - (~dm)(1,2)*tmp41 + (~dm)(1,3)*tmp42 - (~dm)(1,4)*tmp43 + (~dm)(1,5)*tmp44 ) +
          (~dm)(0,2) * ( (~dm)(1,0)*tmp37 - (~dm)(1,1)*tmp41 + (~dm)(1,3)*tmp45 - (~dm)(1,4)*tmp46 + (~dm)(1,5)*tmp47 ) -
          (~dm)(0,3) * ( (~dm)(1,0)*tmp38 - (~dm)(1,1)*tmp42 + (~dm)(1,2)*tmp45 - (~dm)(1,4)*tmp48 + (~dm)(1,5)*tmp49 ) +
          (~dm)(0,4) * ( (~dm)(1,0)*tmp39 - (~dm)(1,1)*tmp43 + (~dm)(1,2)*tmp46 - (~dm)(1,3)*tmp48 + (~dm)(1,5)*tmp50 ) -
          (~dm)(0,5) * ( (~dm)(1,0)*tmp40 - (~dm)(1,1)*tmp44 + (~dm)(1,2)*tmp47 - (~dm)(1,3)*tmp49 + (~dm)(1,4)*tmp50 );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given dense square matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The determinant of the given matrix.
//
// This function computes the determinant of the given dense square matrix via an LU decomposition
// of the matrix.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
typename MT::ElementType detNxN( const DenseMatrix<MT,SO>& dm )
{
   typedef typename MT::ResultType           RT;
   typedef typename MT::ElementType          ET;
   typedef typename RemoveAdaptor<RT>::Type  URT;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( URT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( URT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( URT );

   BLAZE_INTERNAL_ASSERT( isSquare( ~dm ), "Non-square symmetric matrix detected" );

   URT A( ~dm );

   int n   ( boost::numeric_cast<int>( A.rows()      ) );
   int lda ( boost::numeric_cast<int>( A.spacing() ) );
   int info( 0 );

   const UniqueArray<int> ipiv( new int[n] );

   getrf( &n, &n, A.data(), &lda, ipiv.get(), &info );

   if( info > 0 )
      return ET(0);

   ET determinant = ET(1);

   for( int i=0; i<n; ++i ) {
      determinant *= ( ipiv[i] == (i+1) )?( A(i,i) ):( -A(i,i) );
   }

   return determinant;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computation of the determinant of the given dense square matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The determinant of the given matrix.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function computes the determinant of the given dense square matrix. The computation fails
// if the given matrix is not a square matrix. In this case either a compilation error is created
// (if possible) or a \a std::invalid_argument exception is thrown.
//
// \note The computation of the determinant is numerically unreliable since especially for large
// matrices the value can overflow during the computation. Please note that this function does
// not guarantee that it is possible to compute the determinant with the given matrix!
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline typename MT::ElementType det( const DenseMatrix<MT,SO>& dm )
{
   typedef typename MT::ElementType  ET;

   if( !isSquare( ~dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   const size_t N( (~dm).rows() );

   if( IsStrictlyTriangular<MT>::value || N == 0UL )
   {
      return ET(0);
   }
   else if( IsUniTriangular<MT>::value )
   {
      return ET(1);
   }
   else if( N == 1UL )
   {
      return (~dm)(0,0);
   }
   else if( IsTriangular<MT>::value )
   {
      ET determinant( (~dm)(0,0) );

      for( size_t i=1UL; i<N; ++i ) {
         determinant *= (~dm)(i,i);
      }

      return determinant;
   }
   else if( N == 2UL )
   {
      return det2x2( ~dm );
   }
   else if( N == 3UL )
   {
      return det3x3( ~dm );
   }
   else if( N == 4UL )
   {
      return det4x4( ~dm );
   }
   else if( N == 5UL )
   {
      return det5x5( ~dm );
   }
   else if( N == 6UL )
   {
      return det6x6( ~dm );
   }
   else
   {
      return detNxN( ~dm );
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
