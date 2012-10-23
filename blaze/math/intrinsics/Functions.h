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

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/math/intrinsics/DerivedTypes.h>
#include <blaze/system/SSE.h>
#include <blaze/util/Assert.h>


namespace blaze {

//=================================================================================================
//
//  LOAD FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Loads a vector of 'short' values.
// \ingroup intrinsics
//
// \param address The first 'short' value to be loaded.
// \return The loaded vector of 'short' values.
*/
inline sse_short_t load( const short* address )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'unsigned short' values.
// \ingroup intrinsics
//
// \param address The first 'unsigned short' value to be loaded.
// \return The loaded vector of 'unsigned short' values.
*/
inline sse_ushort_t load( const unsigned short* address )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'int' values.
// \ingroup intrinsics
//
// \param address The first 'int' value to be loaded.
// \return The loaded vector of 'int' values.
*/
inline sse_int_t load( const int* address )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'unsigned int' values.
// \ingroup intrinsics
//
// \param address The first 'unsigned int' value to be loaded.
// \return The loaded vector of 'unsigned int' values.
*/
inline sse_uint_t load( const unsigned int* address )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'long' values.
// \ingroup intrinsics
//
// \param address The first 'long' value to be loaded.
// \return The loaded vector of 'long' values.
*/
inline sse_long_t load( const long* address )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'unsigned long' values.
// \ingroup intrinsics
//
// \param address The first 'unsigned long' value to be loaded.
// \return The loaded vector of 'unsigned long' values.
*/
inline sse_ulong_t load( const unsigned long* address )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'float' values.
// \ingroup intrinsics
//
// \param address The first 'float' value to be loaded.
// \return The loaded vector of 'float' values.
*/
inline sse_float_t load( const float* address )
{
#if BLAZE_MIC_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 64UL ), "Invalid alignment detected" );
   return _mm512_load_ps( address );
#elif BLAZE_AVX_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 32UL ), "Invalid alignment detected" );
   return _mm256_load_ps( address );
#elif BLAZE_SSE_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   return _mm_load_ps( address );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'double' values.
// \ingroup intrinsics
//
// \param address The first 'double' value to be loaded.
// \return The loaded vector of 'double' values.
*/
inline sse_double_t load( const double* address )
{
#if BLAZE_MIC_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 64UL ), "Invalid alignment detected" );
   return _mm512_load_pd( address );
#elif BLAZE_AVX_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 32UL ), "Invalid alignment detected" );
   return _mm256_load_pd( address );
#elif BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   return _mm_load_pd( address );
#else
   return *address;
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  SET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Sets all values in the vector to the given 'short' value.
// \ingroup intrinsics
//
// \param value The given 'short' value.
// \return The set vector of 'short' values.
*/
inline sse_short_t set( short value )
{
#if BLAZE_SSE2_MODE
   return _mm_set1_epi16( value );
#else
   return value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets all values in the vector to the given 'unsigned short' value.
// \ingroup intrinsics
//
// \param value The given 'unsigned short' value.
// \return The set vector of 'unsigned short' values.
*/
inline sse_ushort_t set( unsigned short value )
{
#if BLAZE_SSE2_MODE
   return _mm_set1_epi16( value );
#else
   return value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets all values in the vector to the given 'int' value.
// \ingroup intrinsics
//
// \param value The given 'int' value.
// \return The set vector of 'int' values.
*/
inline sse_int_t set( int value )
{
#if BLAZE_SSE2_MODE
   return _mm_set1_epi32( value );
#else
   return value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets all values in the vector to the given 'unsigned int' value.
// \ingroup intrinsics
//
// \param value The given 'unsigned int' value.
// \return The set vector of 'unsigned int' values.
*/
inline sse_uint_t set( unsigned int value )
{
#if BLAZE_SSE2_MODE
   return _mm_set1_epi32( value );
#else
   return value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets all values in the vector to the given 'float' value.
// \ingroup intrinsics
//
// \param value The given 'float' value.
// \return The set vector of 'float' values.
*/
inline sse_float_t set( float value )
{
#if BLAZE_MIC_MODE
   return _mm512_set1_ps( value );
#elif BLAZE_AVX_MODE
   return _mm256_set1_ps( value );
#elif BLAZE_SSE_MODE
   return _mm_set1_ps( value );
#else
   return value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets all values in the vector to the given 'double' value.
// \ingroup intrinsics
//
// \param value The given 'double' value.
// \return The set vector of 'double' values.
*/
inline sse_double_t set( double value )
{
#if BLAZE_MIC_MODE
   return _mm512_set1_pd( value );
#elif BLAZE_AVX_MODE
   return _mm256_set1_pd( value );
#elif BLAZE_SSE2_MODE
   return _mm_set1_pd( value );
#else
   return value;
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  STORE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Aligned store of a vector of 'short' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'short' vector to be stored.
// \return void
*/
inline void store( short* address, sse_short_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_store_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a vector of 'unsigned short' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'unsigned short' vector to be stored.
// \return void
*/
inline void store( unsigned short* address, sse_ushort_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_store_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a vector of 'int' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'int' vector to be stored.
// \return void
*/
inline void store( int* address, sse_int_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_store_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a vector of 'unsigned int' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'unsigned int' vector to be stored.
// \return void
*/
inline void store( unsigned int* address, sse_uint_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_store_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a vector of 'long' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'long' vector to be stored.
// \return void
*/
inline void store( long* address, sse_long_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_store_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a vector of 'unsigned long' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'unsigned long' vector to be stored.
// \return void
*/
inline void store( unsigned long* address, sse_ulong_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_store_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a vector of 'float' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'float' vector to be stored.
// \return void
*/
inline void store( float* address, sse_float_t value )
{
#if BLAZE_MIC_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 64UL ), "Invalid alignment detected" );
   _mm512_store_ps( address, value.value );
#elif BLAZE_AVX_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 32UL ), "Invalid alignment detected" );
   _mm256_store_ps( address, value.value );
#elif BLAZE_SSE_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_store_ps( address, value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a vector of 'double' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'double' vector to be stored.
// \return void
*/
inline void store( double* address, sse_double_t value )
{
#if BLAZE_MIC_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 64UL ), "Invalid alignment detected" );
   _mm512_store_pd( address, value.value );
#elif BLAZE_AVX_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 32UL ), "Invalid alignment detected" );
   _mm256_store_pd( address, value.value );
#elif BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_store_pd( address, value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  STREAM FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'short' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'short' vector to be streamed.
// \return void
*/
inline void stream( short* address, sse_short_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'unsigned short' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'unsigned short' vector to be streamed.
// \return void
*/
inline void stream( unsigned short* address, sse_ushort_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'int' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'int' vector to be streamed.
// \return void
*/
inline void stream( int* address, sse_int_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'unsigned int' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'unsigned int' vector to be streamed.
// \return void
*/
inline void stream( unsigned int* address, sse_uint_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'long' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'long' vector to be streamed.
// \return void
*/
inline void stream( long* address, sse_long_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'unsigned long' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'unsigned long' vector to be streamed.
// \return void
*/
inline void stream( unsigned long* address, sse_ulong_t value )
{
#if BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_stream_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'float' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'float' vector to be streamed.
// \return void
*/
inline void stream( float* address, sse_float_t value )
{
#if BLAZE_MIC_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 64UL ), "Invalid alignment detected" );
   _mm512_stream_ps( address, value.value );
#elif BLAZE_AVX_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 32UL ), "Invalid alignment detected" );
   _mm256_stream_ps( address, value.value );
#elif BLAZE_SSE_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_stream_ps( address, value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a vector of 'double' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'double' vector to be streamed.
// \return void
*/
inline void stream( double* address, sse_double_t value )
{
#if BLAZE_MIC_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 64UL ), "Invalid alignment detected" );
   _mm512_store_pd( address, value.value );
#elif BLAZE_AVX_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 32UL ), "Invalid alignment detected" );
   _mm256_stream_pd( address, value.value );
#elif BLAZE_SSE2_MODE
   BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
   _mm_stream_pd( address, value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************




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
#if BLAZE_SSE2_MODE
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
#if BLAZE_SSE2_MODE
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




//=================================================================================================
//
//  INTRINSIC ADDITION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\fn sse_int8_t operator+( sse_int8_t, sse_int8_t )
// \brief Addition of two vectors of 8-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the addition.
*/
#if BLAZE_SSE2_MODE
inline sse_int8_t operator+( sse_int8_t a, sse_int8_t b )
{
   return _mm_add_epi8( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_int16_t operator+( sse_int16_t, sse_int16_t )
// \brief Addition of two vectors of 16-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the addition.
*/
#if BLAZE_SSE2_MODE
inline sse_int16_t operator+( sse_int16_t a, sse_int16_t b )
{
   return _mm_add_epi16( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_int32_t operator+( sse_int32_t, sse_int32_t )
// \brief Addition of two vectors of 32-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the addition.
*/
#if BLAZE_SSE2_MODE
inline sse_int32_t operator+( sse_int32_t a, sse_int32_t b )
{
   return _mm_add_epi32( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_int64_t operator+( sse_int64_t, sse_int64_t )
// \brief Addition of two vectors of 64-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the addition.
*/
#if BLAZE_SSE2_MODE
inline sse_int64_t operator+( sse_int64_t a, sse_int64_t b )
{
   return _mm_add_epi64( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_float_t operator+( sse_float_t, sse_float_t )
// \brief Addition of two vectors of single precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the addition.
*/
#if BLAZE_MIC_MODE
inline sse_float_t operator+( sse_float_t a, sse_float_t b )
{
   return _mm512_add_ps( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_float_t operator+( sse_float_t a, sse_float_t b )
{
   return _mm256_add_ps( a.value, b.value );
}
#elif BLAZE_SSE_MODE
inline sse_float_t operator+( sse_float_t a, sse_float_t b )
{
   return _mm_add_ps( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_double_t operator+( sse_double_t, sse_double_t )
// \brief Addition of two vectors of double precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the addition.
*/
#if BLAZE_MIC_MODE
inline sse_double_t operator+( sse_double_t a, sse_double_t b )
{
   return _mm512_add_pd( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_double_t operator+( sse_double_t a, sse_double_t b )
{
   return _mm256_add_pd( a.value, b.value );
}
#elif BLAZE_SSE2_MODE
inline sse_double_t operator+( sse_double_t a, sse_double_t b )
{
   return _mm_add_pd( a.value, b.value );
}
#endif
//*************************************************************************************************




//=================================================================================================
//
//  INTRINSIC SUBTRACTION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\fn sse_int8_t operator-( sse_int8_t, sse_int8_t )
// \brief Subtraction of two vectors of 8-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the subtraction.
*/
#if BLAZE_SSE2_MODE
inline sse_int8_t operator-( sse_int8_t a, sse_int8_t b )
{
   return _mm_sub_epi8( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_int16_t operator-( sse_int16_t, sse_int16_t )
// \brief Subtraction of two vectors of 16-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the subtraction.
*/
#if BLAZE_SSE2_MODE
inline sse_int16_t operator-( sse_int16_t a, sse_int16_t b )
{
   return _mm_sub_epi16( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_int32_t operator-( sse_int32_t, sse_int32_t )
// \brief Subtraction of two vectors of 32-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the subtraction.
*/
#if BLAZE_SSE2_MODE
inline sse_int32_t operator-( sse_int32_t a, sse_int32_t b )
{
   return _mm_sub_epi32( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_int64_t operator-( sse_int64_t, sse_int64_t )
// \brief Subtraction of two vectors of 64-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the subtraction.
*/
#if BLAZE_SSE2_MODE
inline sse_int64_t operator-( sse_int64_t a, sse_int64_t b )
{
   return _mm_sub_epi64( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_float_t operator-( sse_float_t, sse_float_t )
// \brief Subtraction of two vectors of single precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the subtraction.
*/
#if BLAZE_MIC_MODE
inline sse_float_t operator-( sse_float_t a, sse_float_t b )
{
   return _mm512_sub_ps( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_float_t operator-( sse_float_t a, sse_float_t b )
{
   return _mm256_sub_ps( a.value, b.value );
}
#elif BLAZE_SSE_MODE
inline sse_float_t operator-( sse_float_t a, sse_float_t b )
{
   return _mm_sub_ps( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_double_t operator-( sse_double_t, sse_double_t )
// \brief Subtraction of two vectors of double precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the subtraction.
*/
#if BLAZE_MIC_MODE
inline sse_double_t operator-( sse_double_t a, sse_double_t b )
{
   return _mm512_sub_pd( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_double_t operator-( sse_double_t a, sse_double_t b )
{
   return _mm256_sub_pd( a.value, b.value );
}
#elif BLAZE_SSE2_MODE
inline sse_double_t operator-( sse_double_t a, sse_double_t b )
{
   return _mm_sub_pd( a.value, b.value );
}
#endif
//*************************************************************************************************




//=================================================================================================
//
//  INTRINSIC MULTIPLICATION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\fn sse_int16_t operator*( sse_int16_t, sse_int16_t )
// \brief Multiplication of two vectors of 16-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the multiplication.
*/
#if BLAZE_SSE2_MODE
inline sse_int16_t operator*( sse_int16_t a, sse_int16_t b )
{
   return _mm_mullo_epi16( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_int32_t operator*( sse_int32_t, sse_int32_t )
// \brief Multiplication of two vectors of 32-bit integral values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the multiplication.
*/
#if BLAZE_SSE4_MODE
inline sse_int32_t operator*( sse_int32_t a, sse_int32_t b )
{
   return _mm_mullo_epi32( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_float_t operator*( sse_float_t, sse_float_t )
// \brief Multiplication of two vectors of single precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the multiplication.
*/
#if BLAZE_MIC_MODE
inline sse_float_t operator*( sse_float_t a, sse_float_t b )
{
   return _mm512_mul_ps( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_float_t operator*( sse_float_t a, sse_float_t b )
{
   return _mm256_mul_ps( a.value, b.value );
}
#elif BLAZE_SSE_MODE
inline sse_float_t operator*( sse_float_t a, sse_float_t b )
{
   return _mm_mul_ps( a.value, b.value );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\fn sse_double_t operator*( sse_double_t, sse_double_t )
// \brief Multiplication of two vectors of double precision floating point values.
// \ingroup intrinsics
//
// \param a The left-hand side operand.
// \param b The right-hand side operand.
// \return The result of the multiplication.
*/
#if BLAZE_MIC_MODE
inline sse_double_t operator*( sse_double_t a, sse_double_t b )
{
   return _mm512_mul_pd( a.value, b.value );
}
#elif BLAZE_AVX_MODE
inline sse_double_t operator*( sse_double_t a, sse_double_t b )
{
   return _mm256_mul_pd( a.value, b.value );
}
#elif BLAZE_SSE2_MODE
inline sse_double_t operator*( sse_double_t a, sse_double_t b )
{
   return _mm_mul_pd( a.value, b.value );
}
#endif
//*************************************************************************************************




//=================================================================================================
//
//  INTRINSIC DIVISION OPERATORS
//
//=================================================================================================

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




//=================================================================================================
//
//  INTRINSIC SUM OPERATION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 16-bit integral intrinsic vector.
// \ingroup intrinsics
//
// \param a The vector to be sumed up.
// \return The sum of all vector elements.
*/
inline int16_t sum( sse_int16_t a )
{
#if BLAZE_SSSE3_MODE
   const sse_int16_t b( _mm_hadd_epi16( a.value, a.value ) );
   const sse_int16_t c( _mm_hadd_epi16( b.value, b.value ) );
   const sse_int16_t d( _mm_hadd_epi16( c.value, c.value ) );
   return d.values[0];
#elif BLAZE_SSE2_MODE
   return a.values[0] + a.values[1] + a.values[2] + a.values[3] +
          a.values[4] + a.values[5] + a.values[6] + a.values[7];
#else
   return a.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the 32-bit integral intrinsic vector.
// \ingroup intrinsics
//
// \param a The vector to be sumed up.
// \return The sum of all vector elements.
*/
inline int32_t sum( sse_int32_t a )
{
#if BLAZE_SSSE3_MODE
   const sse_int32_t b( _mm_hadd_epi32( a.value, a.value ) );
   const sse_int32_t c( _mm_hadd_epi32( b.value, b.value ) );
   return c.values[0];
#elif BLAZE_SSE2_MODE
   return a.values[0] + a.values[1] + a.values[2] + a.values[3];
#else
   return a.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the single precision floating point intrinsic vector.
// \ingroup intrinsics
//
// \param a The vector to be sumed up.
// \return The sum of all vector elements.
*/
inline float sum( sse_float_t a )
{
#if BLAZE_MIC_MODE
   return _mm512_reduce_add_ps( a.value );
#elif BLAZE_AVX_MODE
   const sse_float_t b( _mm256_hadd_ps( a.value, a.value ) );
   const sse_float_t c( _mm256_hadd_ps( b.value, b.value ) );
   const sse_float_t d( _mm256_hadd_ps( c.value, c.value ) );
   return d.values[0];
#elif BLAZE_SSE3_MODE
   const sse_float_t b( _mm_hadd_ps( a.value, a.value ) );
   const sse_float_t c( _mm_hadd_ps( b.value, b.value ) );
   return c.values[0];
#elif BLAZE_SSE_MODE
   return a.values[0] + a.values[1] + a.values[2] + a.values[3];
#else
   return a.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the sum of all elements in the double precision floating point intrinsic vector.
// \ingroup intrinsics
//
// \param a The vector to be sumed up.
// \return The sum of all vector elements.
*/
inline double sum( sse_double_t a )
{
#if BLAZE_MIC_MODE
   return _mm512_reduce_add_pd( a.value );
#elif BLAZE_AVX_MODE
   const sse_double_t b( _mm256_hadd_pd( a.value, a.value ) );
   const sse_double_t c( _mm256_hadd_pd( b.value, b.value ) );
   return c.values[0];
#elif BLAZE_SSE3_MODE
   const sse_double_t b( _mm_hadd_pd( a.value, a.value ) );
   return b.values[0];
#elif BLAZE_SSE2_MODE
   return a.values[0] + a.values[1];
#else
   return a.value;
#endif
}
//*************************************************************************************************

} // namespace blaze

#endif
