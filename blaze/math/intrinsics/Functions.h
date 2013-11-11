//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Functions.h
//  \brief Header file for all intrinsic functions
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

#ifndef _BLAZE_MATH_INTRINSICS_FUNCTIONS_H_
#define _BLAZE_MATH_INTRINSICS_FUNCTIONS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/Addition.h>
#include <blaze/math/intrinsics/Division.h>
#include <blaze/math/intrinsics/Load.h>
#include <blaze/math/intrinsics/Loadu.h>
#include <blaze/math/intrinsics/Multiplication.h>
#include <blaze/math/intrinsics/Reduction.h>
#include <blaze/math/intrinsics/Set.h>
#include <blaze/math/intrinsics/Setzero.h>
#include <blaze/math/intrinsics/Store.h>
#include <blaze/math/intrinsics/Storeu.h>
#include <blaze/math/intrinsics/Stream.h>
#include <blaze/math/intrinsics/Subtraction.h>


namespace blaze {

//=================================================================================================
//
//  INTRINSIC DOT PRODUCT
//
//=================================================================================================

//*************************************************************************************************
/*\brief Dot product of two vectors of single precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the dot product.
*/
// #if BLAZE_SSE4_MODE
// inline float dot( sse_float_t a, sse_float_t b )
// {
//    return _mm_cvtss_f32( _mm_dp_ps( a.value, b.value, 0xF1 ) );
// }
// #elif BLAZE_SSE2_MODE
// inline float dot( sse_float_t a, sse_float_t b )
// {
//    float array[4];
//    store( array, a * b );
//    return array[0] + array[1] + array[2] + array[3];
// }
// #endif
//*************************************************************************************************


//*************************************************************************************************
/*\brief Dot product of two vectors of double precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the dot product.
*/
// #if BLAZE_SSE4_MODE
// inline double dot( sse_double_t a, sse_double_t b )
// {
//    return _mm_cvtsd_f64( _mm_dp_pd( a.value, b.value, 0xF1 ) );
// }
// #elif BLAZE_SSE2_MODE
// inline double dot( sse_double_t a, sse_double_t b )
// {
//    double array[2];
//    store( array, a * b );
//    return array[0] + array[1];
// }
// #endif
//*************************************************************************************************

} // namespace blaze

#endif
