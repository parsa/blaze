//=================================================================================================
/*!
//  \file blaze/math/expressions/TDVecDVecMultExpr.h
//  \brief Header file for the dense vector/dense vector inner product expression
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_EXPRESSIONS_TDVECDVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDVECDVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <stdexcept>
#include <boost/type_traits/remove_reference.hpp>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/Intrinsics.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the dense vector/dense vector scalar multiplication.
// \ingroup dense_vector
*/
template< typename T1    // Type of the left-hand side dense vector
        , typename T2 >  // Type of the right-hand side dense vector
struct TDVecDVecMultExprHelper
{
   //**Type definitions****************************************************************************
   //! Composite type of the left-hand side dense vector expression.
   typedef typename boost::remove_reference< typename T1::CompositeType >::type  CT1;

   //! Composite type of the right-hand side dense vector expression.
   typedef typename boost::remove_reference< typename T2::CompositeType >::type  CT2;
   //**********************************************************************************************

   //**********************************************************************************************
   enum { value = CT1::vectorizable &&
                  CT2::vectorizable &&
                  IsSame< typename CT1::ElementType, typename CT2::ElementType>::value &&
                  IntrinsicTrait< typename CT1::ElementType >::addition &&
                  IntrinsicTrait< typename CT2::ElementType >::multiplication };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default multiplication operator for the scalar product (inner product) of two dense
//        vectors (\f$ s=\vec{a}*\vec{b} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the inner product.
// \param rhs The right-hand side dense vector for the inner product.
// \return The scalar product.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the scalar product (inner product) of two dense vectors:

   \code
   blaze::DynamicVector<double> a, b;
   blaze::double res;
   // ... Resizing and initialization
   res = trans(a) * b;
   \endcode

// The operator returns a scalar value of the higher-order element type of the two involved
// vector element types \a T1::ElementType and \a T2::ElementType. Both vector types \a T1
// and \a T2 as well as the two element types \a T1::ElementType and \a T2::ElementType have
// to be supported by the MultTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1    // Type of the left-hand side dense vector
        , typename T2 >  // Type of the right-hand side dense vector
inline typename DisableIf< TDVecDVecMultExprHelper<T1,T2>,
                           const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type >::Type
   operator*( const DenseVector<T1,true>& lhs, const DenseVector<T2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   typedef typename T1::CompositeType         Lhs;
   typedef typename T2::CompositeType         Rhs;
   typedef typename T1::ElementType           ET1;
   typedef typename T2::ElementType           ET2;
   typedef typename MultTrait<ET1,ET2>::Type  MultType;

   if( (~lhs).size() == 0UL ) return MultType();

   Lhs left ( ~lhs );
   Rhs right( ~rhs );

   MultType sp( left[0UL] * right[0UL] );

   for( size_t i=1UL; i<left.size(); ++i )
      sp += left[i] * right[i];

   return sp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Intrinsic optimized multiplication operator for the scalar product (inner product) of
//        two dense vectors (\f$ s=\vec{a}*\vec{b} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the inner product.
// \param rhs The right-hand side dense vector for the inner product.
// \return The scalar product.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the scalar product (inner product) of two dense vectors:

   \code
   using blaze::columnVector;

   blaze::DynamicVector<double,columnVector> a, b;
   blaze::double res;
   // ... Resizing and initialization
   res = trans(a) * b;
   \endcode

// The operator returns a scalar value of the higher-order element type of the two involved
// vector element types \a T1::ElementType and \a T2::ElementType. Both vector types \a T1
// and \a T2 as well as the two element types \a T1::ElementType and \a T2::ElementType have
// to be supported by the MultTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename T1    // Type of the left-hand side dense vector
        , typename T2 >  // Type of the right-hand side dense vector
inline typename EnableIf< TDVecDVecMultExprHelper<T1,T2>,
                          const typename MultTrait<typename T1::ElementType,typename T2::ElementType>::Type >::Type
   operator*( const DenseVector<T1,true>& lhs, const DenseVector<T2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).size() != (~rhs).size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   typedef typename T1::CompositeType         Lhs;
   typedef typename T2::CompositeType         Rhs;
   typedef typename T1::ElementType           ET1;
   typedef typename T2::ElementType           ET2;
   typedef typename MultTrait<ET1,ET2>::Type  MultType;
   typedef IntrinsicTrait<MultType>           IT;

   if( (~lhs).size() == 0UL ) return MultType();

   Lhs left ( ~lhs );
   Rhs right( ~rhs );

   typename IT::Type xmm1, xmm2, xmm3, xmm4;
   MultType sp( 0 );

   const size_t N  ( left.size() );
   const size_t end( N - N % (IT::size*4UL) );

   for( size_t i=0UL; i<end; i+=IT::size*4UL ) {
      xmm1 = xmm1 + ( left.get(i             ) * right.get(i             ) );
      xmm2 = xmm2 + ( left.get(i+IT::size    ) * right.get(i+IT::size    ) );
      xmm3 = xmm3 + ( left.get(i+IT::size*2UL) * right.get(i+IT::size*2UL) );
      xmm4 = xmm4 + ( left.get(i+IT::size*3UL) * right.get(i+IT::size*3UL) );
   }

   MultType array[IT::size];
   store( array, xmm1 + xmm2 + xmm3 + xmm4 );

   for( size_t i=0UL; i<IT::size; ++i )
      sp += array[i];
   for( size_t i=end; i<N; ++i )
      sp += left[i] * right[i];

   return sp;
}
//*************************************************************************************************

} // namespace blaze

#endif
