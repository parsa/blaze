//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Stream.h
//  \brief Header file for the intrinsic stream functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_STREAM_H_
#define _BLAZE_MATH_INTRINSICS_STREAM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/typetraits/HasSize.h>
#include <blaze/util/typetraits/IsIntegral.h>


namespace blaze {

//=================================================================================================
//
//  INTRINSIC STREAM FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 2-byte integral values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 2-byte integral vector to be streamed.
// \return void
*/
template< typename T >  // Type of the integral value
BLAZE_ALWAYS_INLINE typename EnableIf< And< IsIntegral<T>, HasSize<T,2UL> > >::Type
   stream( T* address, const simd_int16_t& value )
{
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_AVX2_MODE
   _mm256_stream_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 4-byte integral values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 4-byte integral vector to be streamed.
// \return void
*/
template< typename T >  // Type of the integral value
BLAZE_ALWAYS_INLINE typename EnableIf< And< IsIntegral<T>, HasSize<T,4UL> > >::Type
   stream( T* address, const simd_int32_t& value )
{
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_MIC_MODE
   _mm512_store_epi32( address, value.value );
#elif BLAZE_AVX2_MODE
   _mm256_stream_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 8-byte integral values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 8-byte integral vector to be streamed.
// \return void
*/
template< typename T >  // Type of the integral value
BLAZE_ALWAYS_INLINE typename EnableIf< And< IsIntegral<T>, HasSize<T,8UL> > >::Type
   stream( T* address, const simd_int64_t& value )
{
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_MIC_MODE
   _mm512_store_epi64( address, value.value );
#elif BLAZE_AVX2_MODE
   _mm256_stream_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'float' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'float' vector to be streamed.
// \return void
*/
BLAZE_ALWAYS_INLINE void stream( float* address, const simd_float_t& value )
{
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_MIC_MODE
   _mm512_storenr_ps( address, value.value );
#elif BLAZE_AVX_MODE
   _mm256_stream_ps( address, value.value );
#elif BLAZE_SSE_MODE
   _mm_stream_ps( address, value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'double' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'double' vector to be streamed.
// \return void
*/
BLAZE_ALWAYS_INLINE void stream( double* address, const simd_double_t& value )
{
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_MIC_MODE
   _mm512_storenr_pd( address, value.value );
#elif BLAZE_AVX_MODE
   _mm256_stream_pd( address, value.value );
#elif BLAZE_SSE2_MODE
   _mm_stream_pd( address, value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 2-byte integral complex values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 2-byte integral complex vector to be streamed.
// \return void
*/
template< typename T >  // Type of the integral value
BLAZE_ALWAYS_INLINE typename EnableIf< And< IsIntegral<T>, HasSize<T,2UL> > >::Type
   stream( complex<T>* address, const simd_cint16_t& value )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<T> ) == 2UL*sizeof( T ) );
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_AVX2_MODE
   _mm256_stream_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 4-byte integral complex values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 4-byte integral complex vector to be streamed.
// \return void
*/
template< typename T >  // Type of the integral value
BLAZE_ALWAYS_INLINE typename EnableIf< And< IsIntegral<T>, HasSize<T,4UL> > >::Type
   stream( complex<T>* address, const simd_cint32_t& value )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<T> ) == 2UL*sizeof( T ) );
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_MIC_MODE
   _mm512_store_epi32( address, value.value );
#elif BLAZE_AVX2_MODE
   _mm256_stream_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 8-byte integral complex values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 8-byte integral complex vector to be streamed.
// \return void
*/
template< typename T >  // Type of the integral value
BLAZE_ALWAYS_INLINE typename EnableIf< And< IsIntegral<T>, HasSize<T,8UL> > >::Type
   stream( complex<T>* address, const simd_cint64_t& value )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<T> ) == 2UL*sizeof( T ) );
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_MIC_MODE
   _mm512_store_epi64( address, value.value );
#elif BLAZE_AVX2_MODE
   _mm256_stream_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'complex<float>' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'complex<float>' vector to be streamed.
// \return void
*/
BLAZE_ALWAYS_INLINE void stream( complex<float>* address, const simd_cfloat_t& value )
{
   BLAZE_STATIC_ASSERT  ( sizeof( complex<float> ) == 2UL*sizeof( float ) );
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_MIC_MODE
   _mm512_storenr_ps( reinterpret_cast<float*>( address ), value.value );
#elif BLAZE_AVX_MODE
   _mm256_stream_ps( reinterpret_cast<float*>( address ), value.value );
#elif BLAZE_SSE_MODE
   _mm_stream_ps( reinterpret_cast<float*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'complex<double>' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'complex<double>' vector to be streamed.
// \return void
*/
BLAZE_ALWAYS_INLINE void stream( complex<double>* address, const simd_cdouble_t& value )
{
   BLAZE_STATIC_ASSERT  ( sizeof( complex<double> ) == 2UL*sizeof( double ) );
   BLAZE_INTERNAL_ASSERT( checkAlignment( address ), "Invalid alignment detected" );

#if BLAZE_MIC_MODE
   _mm512_storenr_pd( reinterpret_cast<double*>( address ), value.value );
#elif BLAZE_AVX_MODE
   _mm256_stream_pd( reinterpret_cast<double*>( address ), value.value );
#elif BLAZE_SSE2_MODE
   _mm_stream_pd( reinterpret_cast<double*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************

} // namespace blaze

#endif
