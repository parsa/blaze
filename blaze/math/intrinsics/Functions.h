//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Functions.h
//  \brief Header file for all intrinsic functions
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
