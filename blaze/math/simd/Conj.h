//=================================================================================================
/*!
//  \file blaze/math/simd/Conj.h
//  \brief Header file for the SIMD conj functionality
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

#ifndef _BLAZE_MATH_SIMD_CONJ_H_
#define _BLAZE_MATH_SIMD_CONJ_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/simd/BasicTypes.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  INTRINSIC COMPLEX CONJUGATE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Complex conjugate of a vector of 8-bit integral values.
// \ingroup simd
//
// \param a The vector of 8-bit integral values.
// \return The complex conjugate values.
*/
template< typename T >  // Type of the operand
BLAZE_ALWAYS_INLINE const T conj( const simd_i8_t<T>& a ) noexcept
{
   return a;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of 16-bit integral values.
// \ingroup simd
//
// \param a The vector of 16-bit integral values.
// \return The complex conjugate values.
*/
template< typename T >  // Type of the operand
BLAZE_ALWAYS_INLINE const T conj( const simd_i16_t<T>& a ) noexcept
{
   return a;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of 16-bit integral complex values.
// \ingroup simd
//
// \param a The vector of 16-bit integral complex values.
// \return The complex conjugate values.
//
// This operation is only available for SSE2 and AVX2.
*/
BLAZE_ALWAYS_INLINE const simd_cint16_t conj( const simd_cint16_t& a ) noexcept
#if BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi16( a.value, _mm256_set_epi16( -1, 1, -1, 1, -1, 1, -1, 1,
                                                         -1, 1, -1, 1, -1, 1, -1, 1 ) );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mullo_epi16( a.value, _mm_set_epi16( -1, 1, -1, 1, -1, 1, -1, 1 ) );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of 32-bit integral values.
// \ingroup simd
//
// \param a The vector of 32-bit integral values.
// \return The complex conjugate values.
*/
template< typename T >  // Type of the operand
BLAZE_ALWAYS_INLINE const T conj( const simd_i32_t<T>& a ) noexcept
{
   return a;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of 32-bit integral complex values.
// \ingroup simd
//
// \param a The vector of 32-bit integral complex values.
// \return The complex conjugate values.
//
// This operation is only available for SSE4, AVX2, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cint32_t conj( const simd_cint32_t& a ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mullo_epi32( a.value, _mm512_set_epi32( -1, 1, -1, 1, -1, 1, -1, 1,
                                                         -1, 1, -1, 1, -1, 1, -1, 1 ) );
}
#elif BLAZE_AVX2_MODE
{
   return _mm256_mullo_epi32( a.value, _mm256_set_epi32( -1, 1, -1, 1, -1, 1, -1, 1 ) );
}
#elif BLAZE_SSE4_MODE
{
   return _mm_mullo_epi32( a.value, _mm_set_epi32( -1, 1, -1, 1 ) );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of 64-bit integral values.
// \ingroup simd
//
// \param a The vector of 64-bit integral values.
// \return The complex conjugate values.
*/
template< typename T >  // Type of the operand
BLAZE_ALWAYS_INLINE const T conj( const simd_i64_t<T>& a ) noexcept
{
   return a;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of single precision floating point values.
// \ingroup simd
//
// \param a The vector of single precision floating point values.
// \return The complex conjugate values.
*/
BLAZE_ALWAYS_INLINE const simd_float_t conj( const simd_float_t& a ) noexcept
{
   return a;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of single precision complex values.
// \ingroup simd
//
// \param a The vector of single precision complex values.
// \return The complex conjugate values.
//
// This operation is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cfloat_t conj( const simd_cfloat_t& a ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mul_ps( a.value, _mm512_set_ps( -1.0F, 1.0F, -1.0F, 1.0F, -1.0F, 1.0F, -1.0F, 1.0F,
                                                 -1.0F, 1.0F, -1.0F, 1.0F, -1.0F, 1.0F, -1.0F, 1.0F ) );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_mul_ps( a.value, _mm256_set_ps( -1.0F, 1.0F, -1.0F, 1.0F, -1.0F, 1.0F, -1.0F, 1.0F ) );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mul_ps( a.value, _mm_set_ps( -1.0F, 1.0F, -1.0F, 1.0F ) );
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of double precision floating point values.
// \ingroup simd
//
// \param a The vector of double precision floating point values.
// \return The complex conjugate values.
*/
BLAZE_ALWAYS_INLINE const simd_double_t conj( const simd_double_t& a ) noexcept
{
   return a;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Complex conjugate of a vector of double precision complex values.
// \ingroup simd
//
// \param a The vector of double precision complex values.
// \return The complex conjugate values.
//
// This operation is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const simd_cdouble_t conj( const simd_cdouble_t& a ) noexcept
#if BLAZE_MIC_MODE
{
   return _mm512_mul_pd( a.value, _mm512_set_pd( -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0 ) );
}
#elif BLAZE_AVX_MODE
{
   return _mm256_mul_pd( a.value, _mm256_set_pd( -1.0, 1.0, -1.0, 1.0 ) );
}
#elif BLAZE_SSE2_MODE
{
   return _mm_mul_pd( a.value, _mm_set_pd( -1.0, 1.0 ) );
}
#else
= delete;
#endif
//*************************************************************************************************

} // namespace blaze

#endif
