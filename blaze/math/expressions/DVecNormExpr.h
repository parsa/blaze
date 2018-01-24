//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecNormExpr.h
//  \brief Header file for the dense vector norm expression
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECNORMEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECNORMEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/functors/Cbrt.h>
#include <blaze/math/functors/Noop.h>
#include <blaze/math/functors/Pow2.h>
#include <blaze/math/functors/Pow3.h>
#include <blaze/math/functors/Sqrt.h>
#include <blaze/math/functors/UnaryPow.h>
#include <blaze/math/shims/Evaluate.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/system/Optimizations.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FalseType.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/Bool.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Template.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/HasMember.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the dense vector norms.
// \ingroup dense_vector
*/
template< typename VT       // Type of the dense vector
        , typename Power >  // Type of the power operation
struct DVecNormHelper
{
   //**Type definitions****************************************************************************
   //! Composite type of the dense vector expression.
   using CT = RemoveReference_< CompositeType_<VT> >;
   //**********************************************************************************************

   //**SIMD support detection**********************************************************************
   //! Definition of the HasSIMDEnabled type trait.
   BLAZE_CREATE_HAS_DATA_OR_FUNCTION_MEMBER_TYPE_TRAIT( HasSIMDEnabled, simdEnabled );

   //! Definition of the HasLoad type trait.
   BLAZE_CREATE_HAS_DATA_OR_FUNCTION_MEMBER_TYPE_TRAIT( HasLoad, load );

   //! Helper structure for the detection of the SIMD capabilities of the given custom operation.
   struct UseSIMDEnabledFlag {
      enum : bool { value = Power::BLAZE_TEMPLATE simdEnabled< ElementType_<VT> >() };
   };
   //**********************************************************************************************

   //**********************************************************************************************
   enum : bool { value = useOptimizedKernels &&
                         CT::simdEnabled &&
                         If_< HasSIMDEnabled<Power>, UseSIMDEnabledFlag, HasLoad<Power> >::value &&
                         HasSIMDAdd< ElementType_<CT>, ElementType_<CT> >::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default backend implementation of the norm of a dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \return The norm of the given vector.
//
// This function implements the performance optimized norm of a dense vector. Due to the
// explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case vectorization cannot be applied.
*/
template< typename VT      // Type of the dense vector
        , bool TF          // Transpose flag
        , typename Power   // Type of the power operation
        , typename Root >  // Type of the root operation
inline decltype(auto) norm_backend( const DenseVector<VT,TF>& dv, Power power, Root root, FalseType )
{
   using CT = CompositeType_<VT>;
   using ET = ElementType_<VT>;
   using RT = decltype( evaluate( root( std::declval<ET>() ) ) );

   if( (~dv).size() == 0UL ) return RT();

   CT tmp( ~dv );

   const size_t N( tmp.size() );

   ET norm( power( tmp[0UL] ) );
   size_t i( 1UL );

   for( ; (i+4UL) <= N; i+=4UL ) {
      norm += power( tmp[i] ) + power( tmp[i+1UL] ) + power( tmp[i+2UL] ) + power( tmp[i+3UL] );
   }
   for( ; (i+2UL) <= N; i+=2UL ) {
      norm += power( tmp[i] ) + power( tmp[i+1UL] );
   }
   for( ; i<N; ++i ) {
      norm += power( tmp[i] );
   }

   return evaluate( root( norm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized backend implementation of the norm of a dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \return The norm of the given vector.
//
// This function implements the performance optimized norm of a dense vector. Due to the
// explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case vectorization can be applied.
*/
template< typename VT      // Type of the dense vector
        , bool TF          // Transpose flag
        , typename Power   // Type of the power operation
        , typename Root >  // Type of the root operation
inline decltype(auto) norm_backend( const DenseVector<VT,TF>& dv, Power power, Root root, TrueType )
{
   using CT = CompositeType_<VT>;
   using ET = ElementType_<VT>;
   using RT = decltype( evaluate( root( std::declval<ET>() ) ) );

   enum : size_t { SIMDSIZE = SIMDTrait<ET>::size };

   if( (~dv).size() == 0UL ) return RT();

   CT tmp( ~dv );

   const size_t N( tmp.size() );

   constexpr bool remainder( !usePadding || !IsPadded<VT>::value );

   const size_t ipos( ( remainder )?( N & size_t(-SIMDSIZE) ):( N ) );
   BLAZE_INTERNAL_ASSERT( !remainder || ( N - ( N % SIMDSIZE ) ) == ipos, "Invalid end calculation" );

   SIMDTrait_<ET> xmm1, xmm2, xmm3, xmm4;
   size_t i( 0UL );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      xmm1 += power( tmp.load(i             ) );
      xmm2 += power( tmp.load(i+SIMDSIZE    ) );
      xmm3 += power( tmp.load(i+SIMDSIZE*2UL) );
      xmm4 += power( tmp.load(i+SIMDSIZE*3UL) );
   }
   for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL ) {
      xmm1 += power( tmp.load(i         ) );
      xmm2 += power( tmp.load(i+SIMDSIZE) );
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      xmm1 += power( tmp.load(i) );
   }

   ET norm( sum( xmm1 + xmm2 + xmm3 + xmm4 ) );

   for( ; remainder && i<N; ++i ) {
      norm += power( tmp[i] );
   }

   return evaluate( root( norm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computes a custom norm for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \param power The functor for the power operation.
// \param root The functor for the root operation.
// \return The norm of the given dense vector.
//
// This function computes a custom norm of the given dense vector by means of the given functors.
// The following example demonstrates the computation of the L2 norm by means of the blaze::Pow2
// and blaze::Sqrt functors:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double l2 = norm( a, blaze::Pow2(), blaze::Sqrt() );
   \endcode
*/
template< typename VT      // Type of the dense vector
        , bool TF          // Transpose flag
        , typename Power   // Type of the power operation
        , typename Root >  // Type of the root operation
decltype(auto) norm_backend( const DenseVector<VT,TF>& dv, Power power, Root root )
{
   return norm_backend( ~dv, power, root, Bool< DVecNormHelper<VT,Power>::value >() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the L2 norm for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \return The L2 norm of the given dense vector.
//
// This function computes the L2 norm of the given dense vector:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double l2 = norm( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) norm( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return norm_backend( ~dv, Pow2(), Sqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the squared L2 norm for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \return The squared L2 norm of the given dense vector.
//
// This function computes the squared L2 norm of the given dense vector:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double l2 = sqrNorm( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) sqrNorm( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return norm_backend( ~dv, Pow2(), Noop() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the L1 norm for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \return The L1 norm of the given dense vector.
//
// This function computes the L1 norm of the given dense vector:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double l1 = l1Norm( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) l1Norm( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return norm_backend( ~dv, Noop(), Noop() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the L2 norm for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \return The L2 norm of the given dense vector.
//
// This function computes the L2 norm of the given dense vector:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double l2 = l2Norm( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) l2Norm( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return norm_backend( ~dv, Pow2(), Sqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the L3 norm for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \return The L3 norm of the given dense vector.
//
// This function computes the L3 norm of the given dense vector:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double l3 = l3Norm( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) l3Norm( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return norm_backend( ~dv, Pow3(), Cbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the Lp norm for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \param p The norm parameter (p > 0).
// \return The Lp norm of the given dense vector.
//
// This function computes the Lp norm of the given dense vector, where the norm is specified by
// the runtime argument \a p:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double lp = lpNorm( a, 2.3 );
   \endcode

// \note The norm parameter \a p is expected to be larger than 0. This precondition is only checked
// by a user assertion.
*/
template< typename VT    // Type of the dense vector
        , bool TF        // Transpose flag
        , typename ST >  //
decltype(auto) lpNorm( const DenseVector<VT,TF>& dv, ST p )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_USER_ASSERT( !isZero( p ), "Invalid p for Lp norm detected" );

   using ScalarType = MultTrait_< UnderlyingBuiltin_<VT>, decltype( inv( p ) ) >;
   return norm_backend( ~dv, UnaryPow<ScalarType>( p ), UnaryPow<ScalarType>( inv( p ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the maximum norm for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the norm computation.
// \return The maximum norm of the given dense vector.
//
// This function computes the maximum norm of the given dense vector:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double max = maxNorm( a );
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) maxNorm( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return max( abs( ~dv ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
