//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Setzero.h
//  \brief Header file for the intrinisc setzero functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_SETZERO_H_
#define _BLAZE_MATH_INTRINSICS_SETZERO_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  INTRINSIC SETZERO FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Setting an integral intrinsic type with 16 8-bit data values to zero.
// \ingroup intrinsics
//
// \param value The value to be set to zero.
// \return void
*/
BLAZE_ALWAYS_INLINE void setzero( simd_int8_t& value )
{
#if BLAZE_AVX2_MODE
   value.value = _mm256_setzero_si256();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_si128();
#else
   value.value = 0;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an integral intrinsic type with 8 16-bit data values to zero.
// \ingroup intrinsics
//
// \param value The value to be set to zero.
// \return void
*/
BLAZE_ALWAYS_INLINE void setzero( simd_int16_t& value )
{
#if BLAZE_AVX2_MODE
   value.value = _mm256_setzero_si256();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_si128();
#else
   value.value = 0;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an integral intrinsic type with 4 32-bit data values to zero.
// \ingroup intrinsics
//
// \param value The value to be set to zero.
// \return void
*/
BLAZE_ALWAYS_INLINE void setzero( simd_int32_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_epi32();
#elif BLAZE_AVX2_MODE
   value.value = _mm256_setzero_si256();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_si128();
#else
   value.value = 0;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting an integral intrinsic type with 2 64-bit data values to zero.
// \ingroup intrinsics
//
// \param value The value to be set to zero.
// \return void
*/
BLAZE_ALWAYS_INLINE void setzero( simd_int64_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_epi32();
#elif BLAZE_AVX2_MODE
   value.value = _mm256_setzero_si256();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_si128();
#else
   value.value = 0;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting a floating point intrinsic type with 4 32-bit single precision data values to zero.
// \ingroup intrinsics
//
// \param value The value to be set to zero.
// \return void
*/
BLAZE_ALWAYS_INLINE void setzero( simd_float_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_ps();
#elif BLAZE_AVX_MODE
   value.value = _mm256_setzero_ps();
#elif BLAZE_SSE_MODE
   value.value = _mm_setzero_ps();
#else
   value.value = 0.0F;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting a floating point intrinsic type with 4 32-bit double precision data values to zero.
// \ingroup intrinsics
//
// \param value The value to be set to zero.
// \return void
*/
BLAZE_ALWAYS_INLINE void setzero( simd_double_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_pd();
#elif BLAZE_AVX_MODE
   value.value = _mm256_setzero_pd();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_pd();
#else
   value.value = 0.0;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting a floating point intrinsic type with 4 32-bit single precision complex values to zero.
// \ingroup intrinsics
//
// \param value The value to be set to zero.
// \return void
*/
BLAZE_ALWAYS_INLINE void setzero( simd_cfloat_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_ps();
#elif BLAZE_AVX_MODE
   value.value = _mm256_setzero_ps();
#elif BLAZE_SSE_MODE
   value.value = _mm_setzero_ps();
#else
   value.value = 0.0F;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting a floating point intrinsic type with 4 32-bit double precision complex values to zero.
// \ingroup intrinsics
//
// \param value The value to be set to zero.
// \return void
*/
BLAZE_ALWAYS_INLINE void setzero( simd_cdouble_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_pd();
#elif BLAZE_AVX_MODE
   value.value = _mm256_setzero_pd();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_pd();
#else
   value.value = 0.0;
#endif
}
//*************************************************************************************************

} // namespace blaze

#endif
