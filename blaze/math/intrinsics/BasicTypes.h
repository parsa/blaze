//=================================================================================================
/*!
//  \file blaze/math/intrinsics/BasicTypes.h
//  \brief Header file for the basic intrinsic types
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

#ifndef _BLAZE_MATH_INTRINSICS_BASICTYPES_H_
#define _BLAZE_MATH_INTRINSICS_BASICTYPES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/SSE.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  BASIC INTRINSIC TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::sse_int8_t
// \brief Intrinsic type for 8-bit integral data values.
// \ingroup intrinsics
*/
/*! \cond BLAZE_INTERNAL */
#if BLAZE_SSE2_MODE
union sse_int8_t {
   inline sse_int8_t() : value( _mm_setzero_si128() ) {}
   inline sse_int8_t( __m128i v ) : value( v ) {}
   __m128i value;
   int8_t values[16];
};
#else
union sse_int8_t {
   inline sse_int8_t() : value( 0 ) {}
   inline sse_int8_t( int8_t v ) : value( v ) {}
   int8_t value;
   int8_t values[1];
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::sse_int16_t
// \brief Intrinsic type for 16-bit integral data values.
// \ingroup intrinsics
*/
/*! \cond BLAZE_INTERNAL */
#if BLAZE_SSE2_MODE
union sse_int16_t {
   inline sse_int16_t() : value( _mm_setzero_si128() ) {}
   inline sse_int16_t( __m128i v ) : value( v ) {}
   __m128i value;
   int16_t values[8];
};
#else
union sse_int16_t {
   inline sse_int16_t() : value( 0 ) {}
   inline sse_int16_t( int16_t v ) : value( v ) {}
   int16_t value;
   int16_t values[1];
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::sse_int32_t
// \brief Intrinsic type for 32-bit integral data values.
// \ingroup intrinsics
*/
/*! \cond BLAZE_INTERNAL */
#if BLAZE_MIC_MODE
union sse_int32_t {
   inline sse_int32_t() : value( _mm512_setzero_epi32() ) {}
   inline sse_int32_t( __m512i v ) : value( v ) {}
   __m512i value;
   int32_t values[16];
};
#elif BLAZE_SSE2_MODE
union sse_int32_t {
   inline sse_int32_t() : value( _mm_setzero_si128() ) {}
   inline sse_int32_t( __m128i v ) : value( v ) {}
   __m128i value;
   int32_t values[4];
};
#else
union sse_int32_t {
   inline sse_int32_t() : value( 0 ) {}
   inline sse_int32_t( int32_t v ) : value( v ) {}
   int32_t value;
   int32_t values[1];
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::sse_int64_t
// \brief Intrinsic type for 64-bit integral data values.
// \ingroup intrinsics
*/
/*! \cond BLAZE_INTERNAL */
#if BLAZE_MIC_MODE
union sse_int64_t {
   inline sse_int64_t() : value( _mm512_setzero_epi32() ) {}
   inline sse_int64_t( __m512i v ) : value( v ) {}
   __m512i value;
   int64_t values[8];
};
#elif BLAZE_SSE2_MODE
union sse_int64_t {
   inline sse_int64_t() : value( _mm_setzero_si128() ) {}
   inline sse_int64_t( __m128i v ) : value( v ) {}
   __m128i value;
   int64_t values[2];
};
#else
union sse_int64_t {
   inline sse_int64_t() : value( 0 ) {}
   inline sse_int64_t( int64_t v ) : value( v ) {}
   int64_t value;
   int64_t values[1];
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::sse_float_t
// \brief Intrinsic type for 32-bit single precision floating point data values.
// \ingroup intrinsics
*/
/*! \cond BLAZE_INTERNAL */
#if BLAZE_MIC_MODE
union sse_float_t {
   inline sse_float_t() : value( _mm512_setzero_ps() ) {}
   inline sse_float_t( __m512 v ) : value( v ) {}
   __m512 value;
   float values[16];
};
#elif BLAZE_AVX_MODE
union sse_float_t {
   inline sse_float_t() : value( _mm256_setzero_ps() ) {}
   inline sse_float_t( __m256 v ) : value( v ) {}
   __m256 value;
   float values[8];
};
#elif BLAZE_SSE_MODE
union sse_float_t {
   inline sse_float_t() : value( _mm_setzero_ps() ) {}
   inline sse_float_t( __m128 v ) : value( v ) {}
   __m128 value;
   float values[4];
};
#else
union sse_float_t {
   inline sse_float_t() : value( 0.0F ) {}
   inline sse_float_t( float v ) : value( v ) {}
   float value;
   float values[1];
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::sse_double_t
// \brief Intrinsic type for 64-bit double precision floating point data values.
// \ingroup intrinsics
*/
/*! \cond BLAZE_INTERNAL */
#if BLAZE_MIC_MODE
union sse_double_t {
   inline sse_double_t() : value( _mm512_setzero_pd() ) {}
   inline sse_double_t( __m512d v ) : value( v ) {}
   __m512d value;
   double values[8];
};
#elif BLAZE_AVX_MODE
union sse_double_t {
   inline sse_double_t() : value( _mm256_setzero_pd() ) {}
   inline sse_double_t( __m256d v ) : value( v ) {}
   __m256d value;
   double values[4];
};
#elif BLAZE_SSE2_MODE
union sse_double_t {
   inline sse_double_t() : value( _mm_setzero_pd() ) {}
   inline sse_double_t( __m128d v ) : value( v ) {}
   __m128d value;
   double values[2];
};
#else
union sse_double_t {
   inline sse_double_t() : value( 0.0 ) {}
   inline sse_double_t( double v ) : value( v ) {}
   double value;
   double values[1];
};
#endif
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
