//=================================================================================================
/*!
//  \file blaze/math/simd/Reduction.h
//  \brief Header file for the SIMD reduction functionality
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

#ifndef _BLAZE_MATH_SIMD_REDUCTION_H_
#define _BLAZE_MATH_SIMD_REDUCTION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/simd/BasicTypes.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  8-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 8-bit integral SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
template< typename T >  // Type of the SIMD element
BLAZE_ALWAYS_INLINE ValueType_<T> sum( const SIMDi8<T>& a ) noexcept
{
#if BLAZE_AVX512BW_MODE
   return (~a)[ 0] + (~a)[ 1] + (~a)[ 2] + (~a)[ 3] + (~a)[ 4] + (~a)[ 5] + (~a)[ 6] + (~a)[ 7] +
          (~a)[ 8] + (~a)[ 9] + (~a)[10] + (~a)[11] + (~a)[12] + (~a)[13] + (~a)[14] + (~a)[15] +
          (~a)[16] + (~a)[17] + (~a)[18] + (~a)[19] + (~a)[20] + (~a)[21] + (~a)[22] + (~a)[23] +
          (~a)[24] + (~a)[25] + (~a)[26] + (~a)[27] + (~a)[28] + (~a)[29] + (~a)[30] + (~a)[31] +
          (~a)[32] + (~a)[33] + (~a)[34] + (~a)[35] + (~a)[36] + (~a)[37] + (~a)[38] + (~a)[39] +
          (~a)[40] + (~a)[41] + (~a)[42] + (~a)[43] + (~a)[44] + (~a)[45] + (~a)[46] + (~a)[47] +
          (~a)[48] + (~a)[49] + (~a)[50] + (~a)[51] + (~a)[52] + (~a)[53] + (~a)[54] + (~a)[55] +
          (~a)[56] + (~a)[57] + (~a)[58] + (~a)[59] + (~a)[60] + (~a)[61] + (~a)[62] + (~a)[63];
#elif BLAZE_AVX2_MODE
   return (~a)[ 0] + (~a)[ 1] + (~a)[ 2] + (~a)[ 3] + (~a)[ 4] + (~a)[ 5] + (~a)[ 6] + (~a)[ 7] +
          (~a)[ 8] + (~a)[ 9] + (~a)[10] + (~a)[11] + (~a)[12] + (~a)[13] + (~a)[14] + (~a)[15] +
          (~a)[16] + (~a)[17] + (~a)[18] + (~a)[19] + (~a)[20] + (~a)[21] + (~a)[22] + (~a)[23] +
          (~a)[24] + (~a)[25] + (~a)[26] + (~a)[27] + (~a)[28] + (~a)[29] + (~a)[30] + (~a)[31];
#elif BLAZE_SSE2_MODE
   return (~a)[ 0] + (~a)[ 1] + (~a)[ 2] + (~a)[ 3] + (~a)[ 4] + (~a)[ 5] + (~a)[ 6] + (~a)[ 7] +
          (~a)[ 8] + (~a)[ 9] + (~a)[10] + (~a)[11] + (~a)[12] + (~a)[13] + (~a)[14] + (~a)[15];
#else
   return (~a).value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 8-bit integral complex SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
template< typename T >  // Type of the SIMD element
BLAZE_ALWAYS_INLINE const ValueType_<T> sum( const SIMDci8<T>& a ) noexcept
{
#if BLAZE_AVX512BW_MODE
   return complex<int8_t>( (~a)[ 0] + (~a)[ 1] + (~a)[ 2] + (~a)[ 3] + (~a)[ 4] + (~a)[ 5] + (~a)[ 6] + (~a)[ 7] +
                           (~a)[ 8] + (~a)[ 9] + (~a)[10] + (~a)[11] + (~a)[12] + (~a)[13] + (~a)[14] + (~a)[15] +
                           (~a)[16] + (~a)[17] + (~a)[18] + (~a)[19] + (~a)[20] + (~a)[21] + (~a)[22] + (~a)[23] +
                           (~a)[24] + (~a)[25] + (~a)[26] + (~a)[27] + (~a)[28] + (~a)[29] + (~a)[30] + (~a)[31] );
#elif BLAZE_AVX2_MODE
   return complex<int8_t>( (~a)[0] + (~a)[1] + (~a)[ 2] + (~a)[ 3] + (~a)[ 4] + (~a)[ 5] + (~a)[ 6] + (~a)[ 7] +
                           (~a)[8] + (~a)[9] + (~a)[10] + (~a)[11] + (~a)[12] + (~a)[13] + (~a)[14] + (~a)[15] );
#elif BLAZE_SSE2_MODE
   return complex<int8_t>( (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3] + (~a)[4] + (~a)[5] + (~a)[6] + (~a)[7] );
#else
   return (~a).value;
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  16-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 16-bit integral SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
template< typename T >  // Type of the SIMD element
BLAZE_ALWAYS_INLINE ValueType_<T> sum( const SIMDi16<T>& a ) noexcept
{
#if BLAZE_AVX512BW_MODE
   const __m256i low ( _mm512_castsi512_si256( (~a).value ) );
   const __m256i high( _mm512_extracti64x4_epi64( (~a).value, 1 ) );
   const __m256i b   ( _mm256_hadd_epi16( low, high ) );
   const __m256i c   ( _mm256_hadd_epi16( b, b ) );
   const __m256i d   ( _mm256_hadd_epi16( c, c ) );
   const __m256i e   ( _mm256_hadd_epi16( d, d ) );
   const __m128i f   ( _mm_add_epi16( _mm256_extracti128_si256( e, 1 )
                                    , _mm256_castsi256_si128( e ) ) );
   return _mm_extract_epi16( f, 0 );
#elif BLAZE_AVX2_MODE
   const __m256i b( _mm256_hadd_epi16( (~a).value, (~a).value ) );
   const __m256i c( _mm256_hadd_epi16( b, b ) );
   const __m256i d( _mm256_hadd_epi16( c, c ) );
   const __m128i e( _mm_add_epi16( _mm256_extracti128_si256( d, 1 )
                                 , _mm256_castsi256_si128( d ) ) );
   return _mm_extract_epi16( e, 0 );
#elif BLAZE_SSSE3_MODE
   const __m128i b( _mm_hadd_epi16( (~a).value, (~a).value ) );
   const __m128i c( _mm_hadd_epi16( b, b ) );
   const __m128i d( _mm_hadd_epi16( c, c ) );
   return _mm_extract_epi16( d, 0 );
#elif BLAZE_SSE2_MODE
   return (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3] + (~a)[4] + (~a)[5] + (~a)[6] + (~a)[7];
#else
   return (~a).value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 16-bit integral complex SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
template< typename T >  // Type of the SIMD element
BLAZE_ALWAYS_INLINE const ValueType_<T> sum( const SIMDci16<T>& a ) noexcept
{
#if BLAZE_AVX512BW_MODE
   return complex<int16_t>( (~a)[0] + (~a)[1] + (~a)[ 2] + (~a)[ 3] + (~a)[ 4] + (~a)[ 5] + (~a)[ 6] + (~a)[ 7] +
                            (~a)[8] + (~a)[9] + (~a)[10] + (~a)[11] + (~a)[12] + (~a)[13] + (~a)[14] + (~a)[15] );
#elif BLAZE_AVX2_MODE
   return complex<int16_t>( (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3] + (~a)[4] + (~a)[5] + (~a)[6] + (~a)[7] );
#elif BLAZE_SSE2_MODE
   return complex<int16_t>( (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3] );
#else
   return (~a).value;
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 32-bit integral SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
template< typename T >  // Type of the SIMD element
BLAZE_ALWAYS_INLINE ValueType_<T> sum( const SIMDi32<T>& a ) noexcept
{
#if BLAZE_AVX512F_MODE
   const __m256i low ( _mm512_castsi512_si256( (~a).value ) );
   const __m256i high( _mm512_extracti64x4_epi64( (~a).value, 1 ) );
   const __m256i b   ( _mm256_hadd_epi32( low, high ) );
   const __m256i c   ( _mm256_hadd_epi32( b, b ) );
   const __m256i d   ( _mm256_hadd_epi32( c, c ) );
   const __m128i e   ( _mm_add_epi32( _mm256_extracti128_si256( d, 1 )
                                    , _mm256_castsi256_si128( d ) ) );
   return _mm_extract_epi32( e, 0 );
#elif BLAZE_MIC_MODE
   return _mm512_reduce_add_epi32( (~a).value );
#elif BLAZE_AVX2_MODE
   const __m256i b( _mm256_hadd_epi32( (~a).value, (~a).value ) );
   const __m256i c( _mm256_hadd_epi32( b, b ) );
   const __m128i d( _mm_add_epi32( _mm256_extracti128_si256( c, 1 )
                                 , _mm256_castsi256_si128( c ) ) );
   return _mm_extract_epi32( d, 0 );
#elif BLAZE_SSSE3_MODE
   const __m128i b( _mm_hadd_epi32( (~a).value, (~a).value ) );
   const __m128i c( _mm_hadd_epi32( b, b ) );
   return _mm_cvtsi128_si32( c );
#elif BLAZE_SSE2_MODE
   return (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3];
#else
   return (~a).value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 32-bit integral complex SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
template< typename T >  // Type of the SIMD element
BLAZE_ALWAYS_INLINE const ValueType_<T> sum( const SIMDci32<T>& a ) noexcept
{
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
   return complex<int32_t>( (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3] + (~a)[4] + (~a)[5] + (~a)[6] + (~a)[7] );
#elif BLAZE_AVX2_MODE
   return complex<int32_t>( (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3] );
#elif BLAZE_SSE2_MODE
   return complex<int32_t>( (~a)[0] + (~a)[1] );
#else
   return (~a).value;
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  64-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 64-bit integral SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
template< typename T >  // Type of the SIMD element
BLAZE_ALWAYS_INLINE ValueType_<T> sum( const SIMDi64<T>& a ) noexcept
{
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
   return (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3] + (~a)[4] + (~a)[5] + (~a)[6] + (~a)[7];
#elif BLAZE_AVX2_MODE
   return (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3];
#elif BLAZE_SSE2_MODE
   return (~a)[0] + (~a)[1];
#else
   return (~a).value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 64-bit integral complex SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
template< typename T >  // Type of the SIMD element
BLAZE_ALWAYS_INLINE const ValueType_<T> sum( const SIMDci64<T>& a ) noexcept
{
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
   return complex<int64_t>( (~a)[0] + (~a)[1] + (~a)[2] + (~a)[3] );
#elif BLAZE_AVX2_MODE
   return complex<int64_t>( (~a)[0] + (~a)[1] );
#elif BLAZE_SSE2_MODE
   return (~a)[0];
#else
   return (~a).value;
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the sum of all elements in the single precision floating point SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
BLAZE_ALWAYS_INLINE float sum( const SIMDfloat& a ) noexcept
{
#if BLAZE_AVX512F_MODE
   __m512 b( _mm512_shuffle_f32x4( a.value, a.value, 0b11'10'11'10 ) );
   const __m512 c( _mm512_add_ps( b, a.value ) );
   const __m512 d( _mm512_shuffle_f32x4( c, c, 0b01'01'01'01 ) );
   const __m512 e( _mm512_add_ps( d, c ) );
   const __m512 f( _mm512_castsi512_ps( _mm512_shuffle_epi32( _mm512_castps_si512( e ), _MM_PERM_BADC ) ) );
   const __m512 g( _mm512_add_ps( e, f ) );
   const __m512 h( _mm512_castsi512_ps( _mm512_shuffle_epi32( _mm512_castps_si512( g ), _MM_PERM_CDAB ) ) );
   b = _mm512_add_ps( g, h );
   return _mm_cvtss_f32( _mm512_castps512_ps128( b ) );
#elif BLAZE_MIC_MODE
   return _mm512_reduce_add_ps( a.value );
#elif BLAZE_AVX_MODE
   const __m256 b( _mm256_hadd_ps( a.value, a.value ) );
   const __m256 c( _mm256_hadd_ps( b, b ) );
   const __m128 d( _mm_add_ps( _mm256_extractf128_ps( c, 1 ), _mm256_castps256_ps128( c ) ) );
   return _mm_cvtss_f32( d );
#elif BLAZE_SSE3_MODE
   const __m128 b( _mm_hadd_ps( a.value, a.value ) );
   const __m128 c( _mm_hadd_ps( b, b ) );
   return _mm_cvtss_f32( c );
#elif BLAZE_SSE_MODE
   return a[0] + a[1] + a[2] + a[3];
#else
   return a.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the single precision complex SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
BLAZE_ALWAYS_INLINE const complex<float> sum( const SIMDcfloat& a ) noexcept
{
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
   return complex<float>( a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7] );
#elif BLAZE_AVX_MODE
   return complex<float>( a[0] + a[1] + a[2] + a[3] );
#elif BLAZE_SSE_MODE
   return complex<float>( a[0] + a[1] );
#else
   return a.value;
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  64-BIT FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the sum of all elements in the double precision floating point SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
BLAZE_ALWAYS_INLINE double sum( const SIMDdouble& a ) noexcept
{
#if BLAZE_AVX512F_MODE
   __m512d b( _mm512_shuffle_f64x2( a.value, a.value, 0b11'10'11'10 ) );
   const __m512d c( _mm512_add_pd( a.value, b ) );
   const __m512d d( _mm512_permutex_pd( c, 0b01'00'11'10 ) );
   const __m512d e( _mm512_add_pd( c , d ) );
   const __m512d f( _mm512_permutex_pd( e, 0b10'11'00'01 ) );
   b = _mm512_add_pd( e, f );
   return _mm_cvtsd_f64( _mm512_castpd512_pd128( b ) );
#elif BLAZE_MIC_MODE
   return _mm512_reduce_add_pd( a.value );
#elif BLAZE_AVX_MODE
   const __m256d b( _mm256_hadd_pd( a.value, a.value ) );
   const __m128d c( _mm_add_pd( _mm256_extractf128_pd( b, 1 ), _mm256_castpd256_pd128( b ) ) );
   return _mm_cvtsd_f64( c );
#elif BLAZE_SSE3_MODE
   const __m128d b( _mm_hadd_pd( a.value, a.value ) );
   return _mm_cvtsd_f64( b );
#elif BLAZE_SSE2_MODE
   return a[0] + a[1];
#else
   return a.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the double precision complex SIMD vector.
// \ingroup simd
//
// \param a The vector to be summed up.
// \return The sum of all vector elements.
*/
BLAZE_ALWAYS_INLINE const complex<double> sum( const SIMDcdouble& a ) noexcept
{
#if BLAZE_AVX512F_MODE || BLAZE_MIC_MODE
   return complex<double>( a[0] + a[1] + a[2] + a[3] );
#elif BLAZE_AVX_MODE
   return complex<double>( a[0] + a[1] );
#elif BLAZE_SSE2_MODE
   return a[0];
#else
   return a.value;
#endif
}
//*************************************************************************************************

} // namespace blaze

#endif
