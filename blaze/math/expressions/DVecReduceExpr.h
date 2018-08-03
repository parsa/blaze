//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecReduceExpr.h
//  \brief Header file for the dense vector reduce expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECREDUCEEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECREDUCEEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/functors/Add.h>
#include <blaze/math/functors/Mult.h>
#include <blaze/math/SIMD.h>
#include <blaze/system/Compiler.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/Template.h>
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
/*!\brief Auxiliary helper struct for the dense vector reduction operation.
// \ingroup dense_vector
*/
template< typename VT    // Type of the dense vector
        , typename OP >  // Type of the reduction operation
struct DVecReduceExprHelper
{
   //**Type definitions****************************************************************************
   //! Composite type of the dense vector expression.
   using CT = RemoveReference_t< CompositeType_t<VT> >;

   //! Element type of the dense vector expression.
   using ET = ElementType_t<CT>;

   //! Definition of the HasSIMDEnabled type trait.
   BLAZE_CREATE_HAS_DATA_OR_FUNCTION_MEMBER_TYPE_TRAIT( HasSIMDEnabled, simdEnabled );

   //! Definition of the HasLoad type trait.
   BLAZE_CREATE_HAS_DATA_OR_FUNCTION_MEMBER_TYPE_TRAIT( HasLoad, load );
   //**********************************************************************************************

   //**SIMD support detection**********************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper structure for the detection of the SIMD capabilities of the given custom operation.
   struct UseSIMDEnabledFlag {
      static constexpr bool test( bool (*fnc)() ) { return fnc(); }
      static constexpr bool test( bool b ) { return b; }
      static constexpr bool value = test( OP::BLAZE_TEMPLATE simdEnabled<ET,ET> );
   };
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   static constexpr bool value =
      ( CT::simdEnabled &&
        If_t< HasSIMDEnabled_v<OP>, UseSIMDEnabledFlag, HasLoad<OP> >::value );
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
/*!\brief Default backend implementation of the reduction of a dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the reduction computation.
// \param op The reduction operation.
// \return The result of the reduction operation.
//
// This function implements the performance optimized reduction operation for a dense vector. Due
// to the explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case vectorization cannot be applied.
*/
template< typename VT    // Type of the dense vector
        , bool TF        // Transpose flag
        , typename OP >  // Type of the reduction operation
inline auto dvecreduce( const DenseVector<VT,TF>& dv, OP op )
   -> DisableIf_t< DVecReduceExprHelper<VT,OP>::value, ElementType_t<VT> >
{
   using CT = CompositeType_t<VT>;
   using ET = ElementType_t<VT>;

   const size_t N( (~dv).size() );

   if( N == 0UL ) return ET{};
   if( N == 1UL ) return (~dv)[0UL];

   CT tmp( ~dv );

   BLAZE_INTERNAL_ASSERT( tmp.size() == N, "Invalid vector size" );

   ET redux1( tmp[0UL] );
   ET redux2( tmp[1UL] );
   size_t i( 2UL );

   for( ; (i+4UL) <= N; i+=4UL ) {
      redux1 = op( op( redux1, tmp[i    ] ), tmp[i+1UL] );
      redux2 = op( op( redux2, tmp[i+2UL] ), tmp[i+3UL] );
   }
   for( ; (i+2UL) <= N; i+=2UL ) {
      redux1 = op( redux1, tmp[i    ] );
      redux2 = op( redux2, tmp[i+1UL] );
   }
   for( ; i<N; ++i ) {
      redux1 = op( redux1, tmp[i] );
   }

   return op( redux1, redux2 );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized backend implementation of the reduction of a dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the reduction computation.
// \param op The reduction operation.
// \return The result of the reduction operation.
//
// This function implements the performance optimized reduction operation for a dense vector. Due
// to the explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case vectorization can be applied.
*/
template< typename VT    // Type of the dense vector
        , bool TF        // Transpose flag
        , typename OP >  // Type of the reduction operation
inline auto dvecreduce( const DenseVector<VT,TF>& dv, OP op )
   -> EnableIf_t< DVecReduceExprHelper<VT,OP>::value, ElementType_t<VT> >
{
   using CT = CompositeType_t<VT>;
   using ET = ElementType_t<VT>;

   const size_t N( (~dv).size() );

   if( N == 0UL ) return ET{};

   CT tmp( ~dv );

   BLAZE_INTERNAL_ASSERT( tmp.size() == N, "Invalid vector size" );

   constexpr size_t SIMDSIZE = SIMDTrait<ET>::size;

   ET redux{};

   if( N >= SIMDSIZE )
   {
      const size_t ipos( N & size_t(-SIMDSIZE) );
      BLAZE_INTERNAL_ASSERT( ( N - ( N % SIMDSIZE ) ) == ipos, "Invalid end calculation" );

      SIMDTrait_t<ET> xmm1( tmp.load(0UL) );

      if( N >= SIMDSIZE*2UL )
      {
         SIMDTrait_t<ET> xmm2( tmp.load(SIMDSIZE) );
         size_t i( SIMDSIZE*2UL );

         for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL ) {
            xmm1 = op( xmm1, tmp.load(i         ) );
            xmm2 = op( xmm2, tmp.load(i+SIMDSIZE) );
         }
         for( ; i<ipos; i+=SIMDSIZE ) {
            xmm1 = op( xmm1, tmp.load(i) );
         }

         xmm1 = op( xmm1, xmm2 );
      }

      redux = reduce( xmm1, op );

      for( size_t i=ipos; i<N; ++i ) {
         redux = op( redux, tmp[i] );
      }
   }
   else {
      redux = tmp[0UL];
      for( size_t i=1UL; i<N; ++i ) {
         redux = op( redux, tmp[i] );
      }
   }

   return redux;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized backend implementation of the summation of a dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the summation.
// \return The result of the summation.
//
// This function implements the performance optimized summation for a dense vector. Due to
// the explicit application of the SFINAE principle, this function can only be selected by
// the compiler in case vectorization can be applied.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline auto dvecreduce( const DenseVector<VT,TF>& dv, Add /*op*/ )
   -> EnableIf_t< DVecReduceExprHelper<VT,Add>::value, ElementType_t<VT> >
{
   using CT = CompositeType_t<VT>;
   using ET = ElementType_t<VT>;

   const size_t N( (~dv).size() );

   if( N == 0UL ) return ET{};

   CT tmp( ~dv );

   BLAZE_INTERNAL_ASSERT( tmp.size() == N, "Invalid vector size" );

   constexpr bool remainder( !usePadding || !IsPadded_v< RemoveReference_t<CT> > );
   constexpr size_t SIMDSIZE = SIMDTrait<ET>::size;

   ET redux{};

   if( !BLAZE_CLANG_COMPILER && !remainder )
   {
      SIMDTrait_t<ET> xmm1, xmm2;
      size_t i( 0UL );

      for( ; (i+SIMDSIZE) < N; i+=SIMDSIZE*2UL ) {
         xmm1 += tmp.load(i         );
         xmm2 += tmp.load(i+SIMDSIZE);
      }
      if( i < N ) {
         xmm1 += tmp.load(i);
      }

      redux = sum( xmm1 + xmm2 );
   }
   else if( !remainder || N >= SIMDSIZE )
   {
      const size_t ipos( ( remainder )?( N & size_t(-SIMDSIZE) ):( N ) );
      BLAZE_INTERNAL_ASSERT( !remainder || ( N - ( N % SIMDSIZE ) ) == ipos, "Invalid end calculation" );

      SIMDTrait_t<ET> xmm1( tmp.load(0UL) );

      if( remainder ? N >= SIMDSIZE*2UL : N > SIMDSIZE )
      {
         SIMDTrait_t<ET> xmm2( tmp.load(SIMDSIZE) );
         size_t i( SIMDSIZE*2UL );

         for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL ) {
            xmm1 += tmp.load(i         );
            xmm2 += tmp.load(i+SIMDSIZE);
         }
         for( ; i<ipos; i+=SIMDSIZE ) {
            xmm1 += tmp.load(i);
         }

         xmm1 += xmm2;
      }

      redux = sum( xmm1 );

      for( size_t i=ipos; remainder && i<N; ++i ) {
         redux += tmp[i];
      }
   }
   else {
      redux = tmp[0UL];
      for( size_t i=1UL; i<N; ++i ) {
         redux += tmp[i];
      }
   }

   return redux;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Performs a custom reduction operation on the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the reduction computation.
// \param op The reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the given dense vector \a dv by means of the given reduction operation
// \a op:

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   const double sum = reduce( a, Add() );
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified. Thus the
// behavior is non-deterministic if \a op is not associative or not commutative. Also, the
// operation is undefined if the given reduction operation modifies the values.
*/
template< typename VT    // Type of the dense vector
        , bool TF        // Transpose flag
        , typename OP >  // Type of the reduction operation
inline decltype(auto) reduce( const DenseVector<VT,TF>& dv, OP op )
{
   BLAZE_FUNCTION_TRACE;

   return dvecreduce( ~dv, op );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduces the given dense vector by means of addition.
// \ingroup dense_vector
//
// \param dv The given dense vector for the reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the given dense vector \a dv by means of addition:

   \code
   blaze::DynamicVector<int> a{ 1, 2, 3, 4 };
   // ... Resizing and initialization
   const int sum = sum( a );  // Results in 10
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) sum( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return reduce( ~dv, Add() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduces the given dense vector by means of multiplication.
// \ingroup dense_vector
//
// \param dv The given dense vector for the reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the given dense vector \a dv by means of multiplication:

   \code
   blaze::DynamicVector<int> a{ 1, 2, 3, 4 };
   // ... Resizing and initialization
   const int prod = prod( a );  // Results in 24
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline decltype(auto) prod( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   return reduce( ~dv, Mult() );
}
//*************************************************************************************************

} // namespace blaze

#endif
