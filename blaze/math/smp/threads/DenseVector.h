//=================================================================================================
/*!
//  \file blaze/math/smp/threads/DenseVector.h
//  \brief Header file for the C++11/Boost thread-based dense vector SMP implementation
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

#ifndef _BLAZE_MATH_SMP_THREADS_DENSEVECTOR_H_
#define _BLAZE_MATH_SMP_THREADS_DENSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/SMPAssignable.h>
#include <blaze/math/DenseSubvector.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Functions.h>
#include <blaze/math/intrinsics/IntrinsicTrait.h>
#include <blaze/math/smp/ParallelSection.h>
#include <blaze/math/smp/SerialSection.h>
#include <blaze/math/smp/threads/ThreadBackend.h>
#include <blaze/math/SparseSubvector.h>
#include <blaze/math/traits/SubvectorExprTrait.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/system/SMP.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  PLAIN ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP assignment of a dense vector to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP assignment
// of a dense vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF2 >    // Transpose flag of the right-hand side dense vector
void smpAssign_backend( DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   typedef typename VT1::ElementType                         ET1;
   typedef typename VT2::ElementType                         ET2;
   typedef IntrinsicTrait<typename VT1::ElementType>         IT;
   typedef typename SubvectorExprTrait<VT1,aligned>::Type    AlignedTarget;
   typedef typename SubvectorExprTrait<VT1,unaligned>::Type  UnalignedTarget;

   const bool vectorizable( VT1::vectorizable && VT2::vectorizable && IsSame<ET1,ET2>::value );
   const bool lhsAligned  ( (~lhs).isAligned() );
   const bool rhsAligned  ( (~rhs).isAligned() );

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (~lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t equalShare   ( (~lhs).size() / threads + addon );
   const size_t rest         ( equalShare & ( IT::size - 1UL ) );
   const size_t sizePerThread( ( vectorizable && rest )?( equalShare - rest + IT::size ):( equalShare ) );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (~lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (~lhs).size() - index ) );

      if( vectorizable && lhsAligned && rhsAligned ) {
         AlignedTarget target( subvector<aligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleAssign( target, subvector<aligned>( ~rhs, index, size ) );
      }
      else if( vectorizable && lhsAligned ) {
         AlignedTarget target( subvector<aligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleAssign( target, subvector<unaligned>( ~rhs, index, size ) );
      }
      else if( vectorizable && rhsAligned ) {
         UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleAssign( target, subvector<aligned>( ~rhs, index, size ) );
      }
      else {
         UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleAssign( target, subvector<unaligned>( ~rhs, index, size ) );
      }
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP assignment of a sparse vector to a dense
//        vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP assignment
// of a sparse vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF2 >    // Transpose flag of the right-hand side sparse vector
void smpAssign_backend( DenseVector<VT1,TF1>& lhs, const SparseVector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   typedef typename VT1::ElementType                         ET1;
   typedef typename VT2::ElementType                         ET2;
   typedef typename SubvectorExprTrait<VT1,unaligned>::Type  UnalignedTarget;

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (~lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t sizePerThread( (~lhs).size() / threads + addon );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (~lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (~lhs).size() - index ) );
      UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
      TheThreadBackend::scheduleAssign( target, subvector<unaligned>( ~rhs, index, size ) );
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP assignment to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP assignment to a dense vector.
// Due to the explicit application of the SFINAE principle, this function can only be selected by
// the compiler in case both operands are SMP-assignable and the element types of both operands
// are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline typename EnableIf< And< IsDenseVector<VT1>
                             , Or< Not< IsSMPAssignable<VT1> >
                                 , Not< IsSMPAssignable<VT2> > > > >::Type
   smpAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );

   assign( ~lhs, ~rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP assignment to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function performs the C++11/Boost thread-based SMP assignment to a dense vector. Due to
// the explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case both operands are SMP-assignable and the element types of both operands are
// not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline typename EnableIf< And< IsDenseVector<VT1>
                             , IsSMPAssignable<VT1>
                             , IsSMPAssignable<VT2> > >::Type
   smpAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( typename VT2::ElementType );

   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(~rhs).canSMPAssign() ) {
         assign( ~lhs, ~rhs );
      }
      else {
         smpAssign_backend( ~lhs, ~rhs );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDITION ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP addition assignment of a dense vector to a
//        dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP addition
// assignment of a dense vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF2 >    // Transpose flag of the right-hand side dense vector
void smpAddAssign_backend( DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   typedef typename VT1::ElementType                         ET1;
   typedef typename VT2::ElementType                         ET2;
   typedef IntrinsicTrait<typename VT1::ElementType>         IT;
   typedef typename SubvectorExprTrait<VT1,aligned>::Type    AlignedTarget;
   typedef typename SubvectorExprTrait<VT1,unaligned>::Type  UnalignedTarget;

   const bool vectorizable( VT1::vectorizable && VT2::vectorizable && IsSame<ET1,ET2>::value );
   const bool lhsAligned  ( (~lhs).isAligned() );
   const bool rhsAligned  ( (~rhs).isAligned() );

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (~lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t equalShare   ( (~lhs).size() / threads + addon );
   const size_t rest         ( equalShare & ( IT::size - 1UL ) );
   const size_t sizePerThread( ( vectorizable && rest )?( equalShare - rest + IT::size ):( equalShare ) );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (~lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (~lhs).size() - index ) );

      if( vectorizable && lhsAligned && rhsAligned ) {
         AlignedTarget target( subvector<aligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleAddAssign( target, subvector<aligned>( ~rhs, index, size ) );
      }
      else if( vectorizable && lhsAligned ) {
         AlignedTarget target( subvector<aligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleAddAssign( target, subvector<unaligned>( ~rhs, index, size ) );
      }
      else if( vectorizable && rhsAligned ) {
         UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleAddAssign( target, subvector<aligned>( ~rhs, index, size ) );
      }
      else {
         UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleAddAssign( target, subvector<unaligned>( ~rhs, index, size ) );
      }
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP addition assignment of a sparse vector to
//        a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP addition
// assignment of a sparse vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF2 >    // Transpose flag of the right-hand side sparse vector
void smpAddAssign_backend( DenseVector<VT1,TF1>& lhs, const SparseVector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   typedef typename VT1::ElementType                         ET1;
   typedef typename VT2::ElementType                         ET2;
   typedef typename SubvectorExprTrait<VT1,unaligned>::Type  UnalignedTarget;

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (~lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t sizePerThread( (~lhs).size() / threads + addon );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (~lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (~lhs).size() - index ) );
      UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
      TheThreadBackend::scheduleAddAssign( target, subvector<unaligned>( ~rhs, index, size ) );
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP addition assignment to a
//        dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP addition assignment to a
// dense vector. Due to the explicit application of the SFINAE principle, this function can only
// be selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline typename EnableIf< And< IsDenseVector<VT1>
                             , Or< Not< IsSMPAssignable<VT1> >
                                 , Not< IsSMPAssignable<VT2> > > > >::Type
   smpAddAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );

   addAssign( ~lhs, ~rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP addition assignment to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function performs the C++11/Boost thread-based SMP addition assignment to a dense vector.
// Due to the explicit application of the SFINAE principle, this function can only be selected by
// the compiler in case both operands are SMP-assignable and the element types of both operands
// are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline typename EnableIf< And< IsDenseVector<VT1>
                             , IsSMPAssignable<VT1>
                             , IsSMPAssignable<VT2> > >::Type
   smpAddAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( typename VT2::ElementType );

   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(~rhs).canSMPAssign() ) {
         addAssign( ~lhs, ~rhs );
      }
      else {
         smpAddAssign_backend( ~lhs, ~rhs );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRACTION ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP subtraction assignment of a dense vector
//        to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP subtraction
// assignment of a dense vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF2 >    // Transpose flag of the right-hand side dense vector
void smpSubAssign_backend( DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   typedef typename VT1::ElementType                         ET1;
   typedef typename VT2::ElementType                         ET2;
   typedef IntrinsicTrait<typename VT1::ElementType>         IT;
   typedef typename SubvectorExprTrait<VT1,aligned>::Type    AlignedTarget;
   typedef typename SubvectorExprTrait<VT1,unaligned>::Type  UnalignedTarget;

   const bool vectorizable( VT1::vectorizable && VT2::vectorizable && IsSame<ET1,ET2>::value );
   const bool lhsAligned  ( (~lhs).isAligned() );
   const bool rhsAligned  ( (~rhs).isAligned() );

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (~lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t equalShare   ( (~lhs).size() / threads + addon );
   const size_t rest         ( equalShare & ( IT::size - 1UL ) );
   const size_t sizePerThread( ( vectorizable && rest )?( equalShare - rest + IT::size ):( equalShare ) );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (~lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (~lhs).size() - index ) );

      if( vectorizable && lhsAligned && rhsAligned ) {
         AlignedTarget target( subvector<aligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleSubAssign( target, subvector<aligned>( ~rhs, index, size ) );
      }
      else if( vectorizable && lhsAligned ) {
         AlignedTarget target( subvector<aligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleSubAssign( target, subvector<unaligned>( ~rhs, index, size ) );
      }
      else if( vectorizable && rhsAligned ) {
         UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleSubAssign( target, subvector<aligned>( ~rhs, index, size ) );
      }
      else {
         UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleSubAssign( target, subvector<unaligned>( ~rhs, index, size ) );
      }
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP subtraction assignment of a sparse vector
//        to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP subtraction
// assignment of a sparse vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF2 >    // Transpose flag of the right-hand side sparse vector
void smpSubAssign_backend( DenseVector<VT1,TF1>& lhs, const SparseVector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   typedef typename VT1::ElementType                         ET1;
   typedef typename VT2::ElementType                         ET2;
   typedef typename SubvectorExprTrait<VT1,unaligned>::Type  UnalignedTarget;

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (~lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t sizePerThread( (~lhs).size() / threads + addon );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (~lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (~lhs).size() - index ) );
      UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
      TheThreadBackend::scheduleSubAssign( target, subvector<unaligned>( ~rhs, index, size ) );
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP subtraction assignment to a
//        dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP subtraction assignment to a
// dense vector. Due to the explicit application of the SFINAE principle, this function can only
// be selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline typename EnableIf< And< IsDenseVector<VT1>
                             , Or< Not< IsSMPAssignable<VT1> >
                                 , Not< IsSMPAssignable<VT2> > > > >::Type
   smpSubAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );

   subAssign( ~lhs, ~rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP subtraction assignment to a dense
//        vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function performs the C++11/Boost thread-based SMP subtraction assignment to a dense
// vector. Due to the explicit application of the SFINAE principle, this function can only be
// selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline typename EnableIf< And< IsDenseVector<VT1>
                             , IsSMPAssignable<VT1>
                             , IsSMPAssignable<VT2> > >::Type
   smpSubAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( typename VT2::ElementType );

   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(~rhs).canSMPAssign() ) {
         subAssign( ~lhs, ~rhs );
      }
      else {
         smpSubAssign_backend( ~lhs, ~rhs );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTIPLICATION ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP multiplication assignment of a dense vector
//        to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP multiplication
// assignment of a dense vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF2 >    // Transpose flag of the right-hand side dense vector
void smpMultAssign_backend( DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   typedef typename VT1::ElementType                         ET1;
   typedef typename VT2::ElementType                         ET2;
   typedef IntrinsicTrait<typename VT1::ElementType>         IT;
   typedef typename SubvectorExprTrait<VT1,aligned>::Type    AlignedTarget;
   typedef typename SubvectorExprTrait<VT1,unaligned>::Type  UnalignedTarget;

   const bool vectorizable( VT1::vectorizable && VT2::vectorizable && IsSame<ET1,ET2>::value );
   const bool lhsAligned  ( (~lhs).isAligned() );
   const bool rhsAligned  ( (~rhs).isAligned() );

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (~lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t equalShare   ( (~lhs).size() / threads + addon );
   const size_t rest         ( equalShare & ( IT::size - 1UL ) );
   const size_t sizePerThread( ( vectorizable && rest )?( equalShare - rest + IT::size ):( equalShare ) );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (~lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (~lhs).size() - index ) );

      if( vectorizable && lhsAligned && rhsAligned ) {
         AlignedTarget target( subvector<aligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleMultAssign( target, subvector<aligned>( ~rhs, index, size ) );
      }
      else if( vectorizable && lhsAligned ) {
         AlignedTarget target( subvector<aligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleMultAssign( target, subvector<unaligned>( ~rhs, index, size ) );
      }
      else if( vectorizable && rhsAligned ) {
         UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleMultAssign( target, subvector<aligned>( ~rhs, index, size ) );
      }
      else {
         UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
         TheThreadBackend::scheduleMultAssign( target, subvector<unaligned>( ~rhs, index, size ) );
      }
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP multiplication assignment of a sparse
//        vector to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP multiplication
// assignment of a sparse vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF2 >    // Transpose flag of the right-hand side sparse vector
void smpMultAssign_backend( DenseVector<VT1,TF1>& lhs, const SparseVector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   typedef typename VT1::ElementType                         ET1;
   typedef typename VT2::ElementType                         ET2;
   typedef typename SubvectorExprTrait<VT1,unaligned>::Type  UnalignedTarget;

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (~lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t sizePerThread( (~lhs).size() / threads + addon );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (~lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (~lhs).size() - index ) );
      UnalignedTarget target( subvector<unaligned>( ~lhs, index, size ) );
      TheThreadBackend::scheduleMultAssign( target, subvector<unaligned>( ~rhs, index, size ) );
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP multiplication assignment
//        to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP multiplication assignment to
// a dense vector. Due to the explicit application of the SFINAE principle, this function can only
// be selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline typename EnableIf< And< IsDenseVector<VT1>
                             , Or< Not< IsSMPAssignable<VT1> >
                                 , Not< IsSMPAssignable<VT2> > > > >::Type
   smpMultAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );

   multAssign( ~lhs, ~rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP multiplication assignment to a dense
//        vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function performs the C++11/Boost thread-based SMP multiplication assignment to a dense
// vector. Due to the explicit application of the SFINAE principle, this function can only be
// selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline typename EnableIf< And< IsDenseVector<VT1>
                             , IsSMPAssignable<VT1>
                             , IsSMPAssignable<VT2> > >::Type
   smpMultAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( typename VT1::ElementType );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( typename VT2::ElementType );

   BLAZE_INTERNAL_ASSERT( (~lhs).size() == (~rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(~rhs).canSMPAssign() ) {
         multAssign( ~lhs, ~rhs );
      }
      else {
         smpMultAssign_backend( ~lhs, ~rhs );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COMPILE TIME CONSTRAINTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( BLAZE_CPP_THREADS_PARALLEL_MODE || BLAZE_BOOST_THREADS_PARALLEL_MODE );

}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
