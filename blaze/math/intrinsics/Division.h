//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Division.h
//  \brief Header file for the intrinisc division functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_DIVISION_H_
#define _BLAZE_MATH_INTRINSICS_DIVISION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/Vectorization.h>


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
inline sse_int32_t operator/( const sse_int32_t& a, const sse_int32_t& b )
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
inline sse_int64_t operator/( const sse_int64_t& a, const sse_int64_t& b )
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
inline sse_float_t operator/( const sse_float_t& a, const sse_float_t& b )
{
   return _mm512_div_ps( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_float_t operator/( const sse_float_t& a, const sse_float_t& b )
{
   return _mm256_div_ps( a.value, b.value );
}
#elif BLAZE_SSE_MODE
inline sse_float_t operator/( const sse_float_t& a, const sse_float_t& b )
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
inline sse_double_t operator/( const sse_double_t& a, const sse_double_t& b )
{
   return _mm512_div_pd( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_double_t operator/( const sse_double_t& a, const sse_double_t& b )
{
   return _mm256_div_pd( a.value, b.value );
}
#elif BLAZE_SSE2_MODE
inline sse_double_t operator/( const sse_double_t& a, const sse_double_t& b )
{
   return _mm_div_pd( a.value, b.value );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
