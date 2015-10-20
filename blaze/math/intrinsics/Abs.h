//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Abs.h
//  \brief Header file for the intrinisc abs functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_ABS_H_
#define _BLAZE_MATH_INTRINSICS_ABS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  INTRINSIC ABSOLUTE VALUE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\fn simd_int8_t abs( simd_int8_t )
// \brief Absolute value of a vector of 8-bit integral values.
// \ingroup intrinsics
//
// \param a The vector of 8-bit integral values.
// \return The absolute values.
*/
#if BLAZE_AVX2_MODE
BLAZE_ALWAYS_INLINE simd_int8_t abs( const simd_int8_t& a )
{
   return _mm256_abs_epi8( a.value );
}
#elif BLAZE_SSSE3_MODE
BLAZE_ALWAYS_INLINE simd_int8_t abs( const simd_int8_t& a )
{
   return _mm_abs_epi8( a.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn simd_int16_t abs( simd_int16_t )
// \brief Absolute value of a vector of 16-bit integral values.
// \ingroup intrinsics
//
// \param a The vector of 16-bit integral values.
// \return The absolute values.
*/
#if BLAZE_AVX2_MODE
BLAZE_ALWAYS_INLINE simd_int16_t abs( const simd_int16_t& a )
{
   return _mm256_abs_epi16( a.value );
}
#elif BLAZE_SSSE3_MODE
BLAZE_ALWAYS_INLINE simd_int16_t abs( const simd_int16_t& a )
{
   return _mm_abs_epi16( a.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn simd_int32_t abs( simd_int32_t )
// \brief Absolute value of a vector of 32-bit integral values.
// \ingroup intrinsics
//
// \param a The vector of 32-bit integral values.
// \return The absolute values.
*/
#if BLAZE_AVX2_MODE
BLAZE_ALWAYS_INLINE simd_int32_t abs( const simd_int32_t& a )
{
   return _mm256_abs_epi32( a.value );
}
#elif BLAZE_SSSE3_MODE
BLAZE_ALWAYS_INLINE simd_int32_t abs( const simd_int32_t& a )
{
   return _mm_abs_epi32( a.value );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
