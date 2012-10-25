//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Setzero.h
//  \brief Header file for the intrinisc setzero functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_SETZERO_H_
#define _BLAZE_MATH_INTRINSICS_SETZERO_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/SSE.h>


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
inline void setzero( sse_int8_t& value )
{
#if BLAZE_SSE2_MODE
   value.value = _mm_setzero_si128();
#else
   value = 0;
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
inline void setzero( sse_int16_t& value )
{
#if BLAZE_SSE2_MODE
   value.value = _mm_setzero_si128();
#else
   value = 0;
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
inline void setzero( sse_int32_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_epi32();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_si128();
#else
   value = 0;
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
inline void setzero( sse_int64_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_epi32();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_si128();
#else
   value = 0;
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
inline void setzero( sse_float_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_ps();
#elif BLAZE_AVX_MODE
   value.value = _mm256_setzero_ps();
#elif BLAZE_SSE_MODE
   value.value = _mm_setzero_ps();
#else
   value = 0.0F;
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
inline void setzero( sse_double_t& value )
{
#if BLAZE_MIC_MODE
   value.value = _mm512_setzero_pd();
#elif BLAZE_AVX_MODE
   value.value = _mm256_setzero_pd();
#elif BLAZE_SSE2_MODE
   value.value = _mm_setzero_pd();
#else
   value = 0.0;
#endif
}
//*************************************************************************************************

} // namespace blaze

#endif
