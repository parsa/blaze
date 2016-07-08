//=================================================================================================
/*!
//  \file blaze/math/simd/FMA.h
//  \brief Header file for the SIMD fused multiply-add (FMA) functionality
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

#ifndef _BLAZE_MATH_SIMD_FMA_H_
#define _BLAZE_MATH_SIMD_FMA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/simd/Addition.h>
#include <blaze/math/simd/BasicTypes.h>
#include <blaze/math/simd/Multiplication.h>
#include <blaze/math/simd/Subtraction.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  SIMD FMADD OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of 16-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 16-bit integral elements in \a a and \a b and adds the
// packed elements in \a c to the intermediate result. This operation is only available for SSE2
// and AVX2.
*/
template< typename T >  // Type of the operands
BLAZE_ALWAYS_INLINE const T
   fmadd( const SIMDi16<T>& a, const SIMDi16<T>& b, const SIMDi16<T>& c ) noexcept
#if BLAZE_SSE2_MODE || BLAZE_AVX2_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of 16-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 16-bit integral elements in \a a and \a b and adds the
// packed elements in \a c to the intermediate result. This operation is only available for SSE2
// and AVX2.
*/
template< typename T1    // Type of the left-hand side multiplication operand
        , typename T2    // Type of the right-hand side multiplication operand
        , typename T3 >  // Type of the right-hand side addition operand
BLAZE_ALWAYS_INLINE const SIMDuint16
   fmadd( const SIMDi16<T1>& a, const SIMDi16<T2>& b, const SIMDi16<T3>& c ) noexcept
#if BLAZE_SSE2_MODE || BLAZE_AVX2_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of 16-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 16-bit integral complex elements in \a a and \a b and
// adds the packed elements in \a c to the intermediate result. This operation is only available
// for SSE2 and AVX2.
*/
template< typename T >  // Type of the operands
BLAZE_ALWAYS_INLINE const T
   fmadd( const SIMDci16<T>& a, const SIMDci16<T>& b, const SIMDci16<T>& c ) noexcept
#if BLAZE_SSE2_MODE || BLAZE_AVX2_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of 32-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 32-bit integral elements in \a a and \a b and adds the
// packed elements in \a c to the intermediate result. This operation is only available for SSE2
// and AVX2.
*/
template< typename T >  // Type of the operands
BLAZE_ALWAYS_INLINE const T
   fmadd( const SIMDi32<T>& a, const SIMDi32<T>& b, const SIMDi32<T>& c ) noexcept
#if BLAZE_SSE4_MODE || BLAZE_AVX2_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of 32-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 32-bit integral elements in \a a and \a b and adds the
// packed elements in \a c to the intermediate result. This operation is only available for SSE2
// and AVX2.
*/
template< typename T1    // Type of the left-hand side multiplication operand
        , typename T2    // Type of the right-hand side multiplication operand
        , typename T3 >  // Type of the right-hand side addition operand
BLAZE_ALWAYS_INLINE const SIMDuint32
   fmadd( const SIMDi32<T1>& a, const SIMDi32<T2>& b, const SIMDi32<T3>& c ) noexcept
#if BLAZE_SSE4_MODE || BLAZE_AVX2_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of 32-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 32-bit integral complex elements in \a a and \a b and
// adds the packed elements in \a c to the intermediate result. This operation is only available
// for SSE2 and AVX2.
*/
template< typename T >  // Type of the operands
BLAZE_ALWAYS_INLINE const T
   fmadd( const SIMDci32<T>& a, const SIMDci32<T>& b, const SIMDci32<T>& c ) noexcept
#if BLAZE_SSE4_MODE || BLAZE_AVX2_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of single precision SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed single-precision (32-bit) floating-point elements in
// \a a and \a b and adds the packed elements in \a c to the intermediate result. This operation
// is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDfloat
   fmadd( const SIMDfloat& a, const SIMDfloat& b, const SIMDfloat& c ) noexcept
#if BLAZE_FMA_MODE && BLAZE_MIC_MODE
{
   return _mm512_fmadd_ps( a.value, b.value, c.value );
}
#elif BLAZE_FMA_MODE && BLAZE_AVX_MODE
{
   return _mm256_fmadd_ps( a.value, b.value, c.value );
}
#elif BLAZE_FMA_MODE && BLAZE_SSE_MODE
{
   return _mm_fmadd_ps( a.value, b.value, c.value );
}
#elif BLAZE_SSE_MODE || BLAZE_AVX_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of single precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed single-precision (32-bit) complex elements in \a a and
// \a b and adds the packed elements in \a c to the intermediate result. This operation is only
// available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDcfloat
   fmadd( const SIMDcfloat& a, const SIMDcfloat& b, const SIMDcfloat& c ) noexcept
#if BLAZE_SSE_MODE || BLAZE_AVX_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of double precision SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed double-precision (64-bit) floating-point elements in
// \a a and \a b and adds the packed elements in \a c to the intermediate result. This operation
// is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDdouble
   fmadd( const SIMDdouble& a, const SIMDdouble& b, const SIMDdouble& c ) noexcept
#if BLAZE_FMA_MODE && BLAZE_MIC_MODE
{
   return _mm512_fmadd_pd( a.value, b.value, c.value );
}
#elif BLAZE_FMA_MODE && BLAZE_AVX_MODE
{
   return _mm256_fmadd_pd( a.value, b.value, c.value );
}
#elif BLAZE_FMA_MODE && BLAZE_SSE2_MODE
{
   return _mm_fmadd_pd( a.value, b.value, c.value );
}
#elif BLAZE_SSE2_MODE || BLAZE_AVX_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-add of three vectors of double precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD addition operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed double-precision (64-bit) complex elements in \a a and
// \a b and adds the packed elements in \a c to the intermediate result. This operation is only
// available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDcdouble
   fmadd( const SIMDcdouble& a, const SIMDcdouble& b, const SIMDcdouble& c ) noexcept
#if BLAZE_SSE2_MODE || BLAZE_AVX_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) + c;
}
#else
= delete;
#endif
//*************************************************************************************************




//=================================================================================================
//
//  SIMD FMSUB OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of 16-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 16-bit integral elements in \a a and \a b and subtracts
// the packed elements in \a c from the intermediate result. This operation is only available for
// SSE2 and AVX2.
*/
template< typename T >  // Type of the operands
BLAZE_ALWAYS_INLINE const T
   fmsub( const SIMDi16<T>& a, const SIMDi16<T>& b, const SIMDi16<T>& c ) noexcept
#if BLAZE_SSE2_MODE || BLAZE_AVX2_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of 16-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 16-bit integral elements in \a a and \a b and subtracts
// the packed elements in \a c from the intermediate result. This operation is only available for
// SSE2 and AVX2.
*/
template< typename T1    // Type of the left-hand side multiplication operand
        , typename T2    // Type of the right-hand side multiplication operand
        , typename T3 >  // Type of the right-hand side subtraction operand
BLAZE_ALWAYS_INLINE const SIMDuint16
   fmsub( const SIMDi16<T1>& a, const SIMDi16<T2>& b, const SIMDi16<T3>& c ) noexcept
#if BLAZE_SSE2_MODE || BLAZE_AVX2_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of 16-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 16-bit integral complex elements in \a a and \a b and
// subtracts the packed elements in \a c from the intermediate result. This operation is only
// available for SSE2 and AVX2.
*/
template< typename T >  // Type of the operands
BLAZE_ALWAYS_INLINE const T
   fmsub( const SIMDci16<T>& a, const SIMDci16<T>& b, const SIMDci16<T>& c ) noexcept
#if BLAZE_SSE2_MODE || BLAZE_AVX2_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of 32-bit integral SIMD values of the same type.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 32-bit integral elements in \a a and \a b and subtracts
// the packed elements in \a c from the intermediate result. This operation is only available for
// SSE2 and AVX2.
*/
template< typename T >  // Type of the operands
BLAZE_ALWAYS_INLINE const T
   fmsub( const SIMDi32<T>& a, const SIMDi32<T>& b, const SIMDi32<T>& c ) noexcept
#if BLAZE_SSE4_MODE || BLAZE_AVX2_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of 32-bit integral SIMD values of different type.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 32-bit integral elements in \a a and \a b and subtracts
// the packed elements in \a c from the intermediate result. This operation is only available for
// SSE2 and AVX2.
*/
template< typename T1    // Type of the left-hand side multiplication operand
        , typename T2    // Type of the right-hand side multiplication operand
        , typename T3 >  // Type of the right-hand side subtraction operand
BLAZE_ALWAYS_INLINE const SIMDuint32
   fmsub( const SIMDi32<T1>& a, const SIMDi32<T2>& b, const SIMDi32<T3>& c ) noexcept
#if BLAZE_SSE4_MODE || BLAZE_AVX2_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of 32-bit integral complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed 32-bit integral complex elements in \a a and \a b and
// subtracts the packed elements in \a c from the intermediate result. This operation is only
// available for SSE2 and AVX2.
*/
template< typename T >  // Type of the operands
BLAZE_ALWAYS_INLINE const T
   fmsub( const SIMDci32<T>& a, const SIMDci32<T>& b, const SIMDci32<T>& c ) noexcept
#if BLAZE_SSE4_MODE || BLAZE_AVX2_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of single precision SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed single-precision (32-bit) floating-point elements in
// \a a and \a b and subtracts the packed elements in \a c from the intermediate result. This
// operation is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDfloat
   fmsub( const SIMDfloat& a, const SIMDfloat& b, const SIMDfloat& c ) noexcept
#if BLAZE_FMA_MODE && BLAZE_MIC_MODE
{
   return _mm512_fmsub_ps( a.value, b.value, c.value );
}
#elif BLAZE_FMA_MODE && BLAZE_AVX_MODE
{
   return _mm256_fmsub_ps( a.value, b.value, c.value );
}
#elif BLAZE_FMA_MODE && BLAZE_SSE_MODE
{
   return _mm_fmsub_ps( a.value, b.value, c.value );
}
#elif BLAZE_SSE_MODE || BLAZE_AVX_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of single precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed single-precision (32-bit) complex elements in \a a and
// \a b and subtracts the packed elements in \a c from the intermediate result. This operation
// is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDcfloat
   fmsub( const SIMDcfloat& a, const SIMDcfloat& b, const SIMDcfloat& c ) noexcept
#if BLAZE_SSE_MODE || BLAZE_AVX_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of double precision SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed double-precision (64-bit) floating-point elements in
// \a a and \a b and subtracts the packed elements in \a c from the intermediate result. This
// operation is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDdouble
   fmsub( const SIMDdouble& a, const SIMDdouble& b, const SIMDdouble& c ) noexcept
#if BLAZE_FMA_MODE && BLAZE_MIC_MODE
{
   return _mm512_fmsub_pd( a.value, b.value, c.value );
}
#elif BLAZE_FMA_MODE && BLAZE_AVX_MODE
{
   return _mm256_fmsub_pd( a.value, b.value, c.value );
}
#elif BLAZE_FMA_MODE && BLAZE_SSE2_MODE
{
   return _mm_fmsub_pd( a.value, b.value, c.value );
}
#elif BLAZE_SSE2_MODE || BLAZE_AVX_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Fused multiply-subtract of three vectors of double precision complex SIMD values.
// \ingroup simd
//
// \param a The left-hand side SIMD multiplication operand.
// \param b The right-hand side SIMD multiplication operand.
// \param c The right-hand side SIMD subtraction operand.
// \return The result of the FMA operation.
//
// This operation multiplies the packed double-precision (64-bit) complex elements in \a a and
// \a b and subtracts the packed elements in \a c from the intermediate result. This operation
// is only available for SSE2, AVX, and AVX-512.
*/
BLAZE_ALWAYS_INLINE const SIMDcdouble
   fmsub( const SIMDcdouble& a, const SIMDcdouble& b, const SIMDcdouble& c ) noexcept
#if BLAZE_SSE2_MODE || BLAZE_AVX_MODE || BLAZE_MIC_MODE
{
   return ( a * b ) - c;
}
#else
= delete;
#endif
//*************************************************************************************************

} // namespace blaze

#endif
