//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Division.h
//  \brief Header file for the intrinisc division functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_DIVISION_H_
#define _BLAZE_MATH_INTRINSICS_DIVISION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/SSE.h>


namespace blaze {

//=================================================================================================
//
//  INTRINSIC DIVISION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\fn sse_float_t operator/( sse_int32_t, sse_int32_t )
// \brief Division of two vectors of 32-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the division.
*/
#if BLAZE_MIC_MODE
inline sse_int32_t operator/( sse_int32_t a, sse_int32_t b )
{
   return _mm512_div_epi32( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_float_t operator/( sse_int64_t, sse_int64_t )
// \brief Division of two vectors of 64-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the division.
*/
#if BLAZE_MIC_MODE
inline sse_int64_t operator/( sse_int64_t a, sse_int64_t b )
{
   return _mm512_div_epi64( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_float_t operator/( sse_float_t, sse_float_t )
// \brief Division of two vectors of single precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the division.
*/
#if BLAZE_MIC_MODE
inline sse_float_t operator/( sse_float_t a, sse_float_t b )
{
   return _mm512_div_ps( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_float_t operator/( sse_float_t a, sse_float_t b )
{
   return _mm256_div_ps( a.value, b.value );
}
#elif BLAZE_SSE_MODE
inline sse_float_t operator/( sse_float_t a, sse_float_t b )
{
   return _mm_div_ps( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_double_t operator/( sse_double_t, sse_double_t )
// \brief Division of two vectors of double precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the division.
*/
#if BLAZE_MIC_MODE
inline sse_double_t operator/( sse_double_t a, sse_double_t b )
{
   return _mm512_div_pd( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_double_t operator/( sse_double_t a, sse_double_t b )
{
   return _mm256_div_pd( a.value, b.value );
}
#elif BLAZE_SSE2_MODE
inline sse_double_t operator/( sse_double_t a, sse_double_t b )
{
   return _mm_div_pd( a.value, b.value );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
