//=================================================================================================
/*!
//  \file blaze/math/simd/BasicTypes.h
//  \brief Header file for the basic SIMD types
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

#ifndef _BLAZE_MATH_SIMD_BASICTYPES_H_
#define _BLAZE_MATH_SIMD_BASICTYPES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsSIMDType.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Vectorization.h>
#include <blaze/util/Complex.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  SIMD BASE TYPES
//
//=================================================================================================

//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_t {
   BLAZE_ALWAYS_INLINE T&       operator~() noexcept       { return *static_cast<T*>( this ); }
   BLAZE_ALWAYS_INLINE const T& operator~() const noexcept { return *static_cast<const T*>( this ); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all 8-bit integral SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_i8_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all 8-bit integral complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_ci8_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all 16-bit integral SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_i16_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all 16-bit integral complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_ci16_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all 32-bit integral SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_i32_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all 32-bit integral complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_ci32_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all 64-bit integral SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_i64_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all 64-bit integral complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_ci64_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all single precision floating point SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_f32_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all single precision floating point complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_cf32_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all double precision floating point SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_f64_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::simd_t
// \brief Base class for all double precision floating point complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct simd_cf64_t : public simd_t< T >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  8-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_int8_t
// \brief SIMD type for 8-bit signed integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_int8_t : public simd_i8_t< simd_int8_t >
{
   using Type = int8_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_int8_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_int8_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_int8_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_int8_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 32 8-bit signed integral data values
   enum : size_t { size = 32UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_int8_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_int8_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 16 8-bit signed integral data values
   enum : size_t { size = 16UL };
#else
   BLAZE_ALWAYS_INLINE simd_int8_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_int8_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_int8_t( const simd_i8_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_int8_t& operator=( const simd_i8_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::simd_uint8_t
// \brief SIMD type for 8-bit unsigned integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_uint8_t : public simd_i8_t< simd_uint8_t >
{
   using Type = uint8_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_uint8_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_uint8_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_uint8_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_uint8_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 32 8-bit unsigned integral data values
   enum : size_t { size = 32UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_uint8_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_uint8_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 16 8-bit unsigned integral data values
   enum : size_t { size = 16UL };
#else
   BLAZE_ALWAYS_INLINE simd_uint8_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_uint8_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_uint8_t( const simd_i8_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_uint8_t& operator=( const simd_i8_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  8-BIT INTEGRAL COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_cint8_t
// \brief SIMD type for 8-bit signed integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cint8_t : public simd_ci8_t< simd_cint8_t >
{
   using Type = complex<int8_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cint8_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cint8_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_cint8_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_cint8_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 16 8-bit signed integral complex values
   enum : size_t { size = 16UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cint8_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_cint8_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 8 8-bit signed integral complex values
   enum : size_t { size = 8UL };
#else
   BLAZE_ALWAYS_INLINE simd_cint8_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cint8_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cint8_t( const simd_ci8_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cint8_t& operator=( const simd_ci8_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::simd_cuint8_t
// \brief SIMD type for 8-bit unsigned integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cuint8_t : public simd_ci8_t< simd_cuint8_t >
{
   using Type = complex<uint8_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cuint8_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cuint8_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_cuint8_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint8_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 16 8-bit unsigned integral complex values
   enum : size_t { size = 16UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cuint8_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint8_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 8 8-bit unsigend integral complex values
   enum : size_t { size = 8UL };
#else
   BLAZE_ALWAYS_INLINE simd_cuint8_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cuint8_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cuint8_t( const simd_ci8_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cuint8_t& operator=( const simd_ci8_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  16-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_int16_t
// \brief SIMD type for 16-bit signed integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_int16_t : public simd_i16_t< simd_int16_t >
{
   using Type = int16_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_int16_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_int16_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_int16_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_int16_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 16 16-bit signed integral data values
   enum : size_t { size = 16UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_int16_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_int16_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 8 16-bit signed integral data values
   enum : size_t { size = 8UL };
#else
   BLAZE_ALWAYS_INLINE simd_int16_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_int16_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_int16_t( const simd_i16_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_int16_t& operator=( const simd_i16_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::simd_uint16_t
// \brief SIMD type for 16-bit unsigned integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_uint16_t : public simd_i16_t< simd_uint16_t >
{
   using Type = uint16_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_uint16_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_uint16_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_uint16_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_uint16_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 16 16-bit unsigned integral data values
   enum : size_t { size = 16UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_uint16_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_uint16_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 8 16-bit unsigned integral data values
   enum : size_t { size = 8UL };
#else
   BLAZE_ALWAYS_INLINE simd_uint16_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_uint16_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_uint16_t( const simd_i16_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_uint16_t& operator=( const simd_i16_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  16-BIT INTEGRAL COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_cint16_t
// \brief SIMD type for 16-bit signed integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cint16_t : public simd_ci16_t< simd_cint16_t >
{
   using Type = complex<int16_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cint16_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cint16_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_cint16_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_cint16_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 8 16-bit signed integral complex values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cint16_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_cint16_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 4 16-bit signed integral complex values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE simd_cint16_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cint16_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cint16_t( const simd_ci16_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cint16_t& operator=( const simd_ci16_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::simd_cuint16_t
// \brief SIMD type for 16-bit unsigned integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cuint16_t : public simd_ci16_t< simd_cuint16_t >
{
   using Type = complex<uint16_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cuint16_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cuint16_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_cuint16_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint16_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 8 16-bit unsigned integral complex values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cuint16_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint16_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 4 16-bit unsigned integral complex values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE simd_cuint16_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cuint16_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cuint16_t( const simd_ci16_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cuint16_t& operator=( const simd_ci16_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_int32_t
// \brief SIMD type for 32-bit signed integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_int32_t : public simd_i32_t< simd_int32_t >
{
   using Type = int32_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_int32_t() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE simd_int32_t( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512i value;  // Contains 16 32-bit signed integral data values
   enum : size_t { size = 16UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_int32_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_int32_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 8 32-bit signed integral data values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_int32_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_int32_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 4 32-bit signed integral data values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE simd_int32_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_int32_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_int32_t( const simd_i32_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_int32_t& operator=( const simd_i32_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::simd_uint32_t
// \brief SIMD type for 32-bit unsigned integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_uint32_t : public simd_i32_t< simd_uint32_t >
{
   using Type = uint32_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_uint32_t() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE simd_uint32_t( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512i value;  // Contains 16 32-bit unsigned integral data values
   enum : size_t { size = 16UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_uint32_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_uint32_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 8 32-bit unsigned integral data values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_uint32_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_uint32_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 4 32-bit unsigned integral data values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE simd_uint32_t() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE simd_uint32_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_uint32_t( const simd_i32_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_uint32_t& operator=( const simd_i32_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT INTEGRAL COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_cint32_t
// \brief SIMD type for 32-bit signed integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cint32_t : public simd_ci32_t< simd_cint32_t >
{
   using Type = complex<int32_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cint32_t() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE simd_cint32_t( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512i value;  // Contains 8 32-bit signed integral complex values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_cint32_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_cint32_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 4 32-bit signed integral complex values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cint32_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_cint32_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 2 32-bit signed integral complex values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE simd_cint32_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cint32_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cint32_t( const simd_ci32_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cint32_t& operator=( const simd_ci32_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::simd_cuint32_t
// \brief SIMD type for 32-bit unsigned integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cuint32_t : public simd_ci32_t< simd_cuint32_t >
{
   using Type = complex<uint32_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cuint32_t() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint32_t( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512i value;  // Contains 8 32-bit unsigned integral complex values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_cuint32_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint32_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 4 32-bit unsigned integral complex values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cuint32_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint32_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 2 32-bit unsigned integral complex values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE simd_cuint32_t() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE simd_cuint32_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cuint32_t( const simd_ci32_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cuint32_t& operator=( const simd_ci32_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  64-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_int64_t
// \brief SIMD type for 64-bit integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_int64_t : public simd_i64_t< simd_int64_t >
{
   using Type = int64_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_int64_t() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE simd_int64_t( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512i value;  // Contains 8 64-bit signed integral data values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_int64_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_int64_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 4 64-bit signed integral data values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_int64_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_int64_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 2 64-bit signed integral data values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE simd_int64_t() noexcept : value( 0L ) {}
   BLAZE_ALWAYS_INLINE simd_int64_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_int64_t( const simd_i64_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_int64_t& operator=( const simd_i64_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::simd_uint64_t
// \brief SIMD type for 64-bit unsigned integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_uint64_t : public simd_i64_t< simd_uint64_t >
{
   using Type = uint64_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_uint64_t() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE simd_uint64_t( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512i value;  // Contains 8 64-bit unsigned integral data values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_uint64_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_uint64_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 4 64-bit unsigned integral data values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_uint64_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_uint64_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 2 64-bit unsigned integral data values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE simd_uint64_t() noexcept : value( 0L ) {}
   BLAZE_ALWAYS_INLINE simd_uint64_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_uint64_t( const simd_i64_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_uint64_t& operator=( const simd_i64_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  64-BIT INTEGRAL COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_cint64_t
// \brief SIMD type for 64-bit signed integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cint64_t : public simd_ci64_t< simd_cint64_t >
{
   using Type = complex<int64_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cint64_t() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE simd_cint64_t( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512i value;  // Contains 4 64-bit signed integral complex values
   enum : size_t { size = 4UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_cint64_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_cint64_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 2 64-bit signed integral complex values
   enum : size_t { size = 2UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cint64_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_cint64_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 1 64-bit signed integral complex values
   enum : size_t { size = 1UL };
#else
   BLAZE_ALWAYS_INLINE simd_cint64_t() noexcept : value( 0L, 0L ) {}
   BLAZE_ALWAYS_INLINE simd_cint64_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cint64_t( const simd_ci64_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cint64_t& operator=( const simd_ci64_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::simd_cuint64_t
// \brief SIMD type for 64-bit unsigned integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cuint64_t : public simd_ci64_t< simd_cuint64_t >
{
   using Type = complex<uint64_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cuint64_t() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint64_t( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512i value;  // Contains 4 64-bit unsigned integral complex values
   enum : size_t { size = 4UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE simd_cuint64_t() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint64_t( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256i value;  // Contains 2 64-bit unsigned integral complex values
   enum : size_t { size = 2UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cuint64_t() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE simd_cuint64_t( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128i value;  // Contains 1 64-bit unsigned integral complex values
   enum : size_t { size = 1UL };
#else
   BLAZE_ALWAYS_INLINE simd_cuint64_t() noexcept : value( 0L, 0L ) {}
   BLAZE_ALWAYS_INLINE simd_cuint64_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cuint64_t( const simd_ci64_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cuint64_t& operator=( const simd_ci64_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SINGLE PRECISION FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_float_t
// \brief SIMD type for 32-bit single precision floating point data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_float_t : public simd_f32_t< simd_float_t >
{
   using Type = float;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_float_t() noexcept : value( _mm512_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE simd_float_t( __m512 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512 value;  // Contains 16 32-bit single precision floating point values
   enum : size_t { size = 16UL };
#elif BLAZE_AVX_MODE
   BLAZE_ALWAYS_INLINE simd_float_t() noexcept : value( _mm256_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE simd_float_t( __m256 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256 value;  // Contains 8 32-bit single precision floating point values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE_MODE
   BLAZE_ALWAYS_INLINE simd_float_t() noexcept : value( _mm_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE simd_float_t( __m128 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128 value;  // Contains 4 32-bit single precision floating point values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE simd_float_t() noexcept : value( 0.0F ) {}
   BLAZE_ALWAYS_INLINE simd_float_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_float_t( const simd_f32_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_float_t& operator=( const simd_f32_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SINGLE PRECISION FLOATING POINT COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_cfloat_t
// \brief SIMD type for 32-bit single precision complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cfloat_t : public simd_cf32_t< simd_cfloat_t >
{
   using Type = complex<float>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cfloat_t() noexcept : value( _mm512_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE simd_cfloat_t( __m512 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512 value;  // Contains 8 32-bit single precision complex values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX_MODE
   BLAZE_ALWAYS_INLINE simd_cfloat_t() noexcept : value( _mm256_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE simd_cfloat_t( __m256 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256 value;  // Contains 4 32-bit single precision complex values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE_MODE
   BLAZE_ALWAYS_INLINE simd_cfloat_t() noexcept : value( _mm_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE simd_cfloat_t( __m128 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128 value;  // Contains 2 32-bit single precision complex values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE simd_cfloat_t() noexcept : value( 0.0F, 0.0F ) {}
   BLAZE_ALWAYS_INLINE simd_cfloat_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cfloat_t( const simd_cf32_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cfloat_t& operator=( const simd_cf32_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DOUBLE PRECISION FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_double_t
// \brief SIMD type for 64-bit double precision floating point data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_double_t : public simd_f64_t< simd_double_t >
{
   using Type = double;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_double_t() noexcept : value( _mm512_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE simd_double_t( __m512d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512d value;  // Contains 8 64-bit double precision floating point values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX_MODE
   BLAZE_ALWAYS_INLINE simd_double_t() noexcept : value( _mm256_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE simd_double_t( __m256d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256d value;  // Contains 4 64-bit double precision floating point values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_double_t() noexcept : value( _mm_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE simd_double_t( __m128d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128d value;  // Contains 2 64-bit double precision floating point values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE simd_double_t() noexcept : value( 0.0 ) {}
   BLAZE_ALWAYS_INLINE simd_double_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_double_t( const simd_f64_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_double_t& operator=( const simd_f64_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DOUBLE PRECISION FLOATING POINT COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::simd_cdouble_t
// \brief SIMD type for 64-bit double precision complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct simd_cdouble_t : public simd_cf64_t< simd_cdouble_t >
{
   using Type = complex<double>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE simd_cdouble_t() noexcept : value( _mm512_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE simd_cdouble_t( __m512d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m512d value;  // Contains 4 64-bit double precision complex value
   enum : size_t { size = 4UL };
#elif BLAZE_AVX_MODE
   BLAZE_ALWAYS_INLINE simd_cdouble_t() noexcept : value( _mm256_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE simd_cdouble_t( __m256d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m256d value;  // Contains 2 64-bit double precision complex value
   enum : size_t { size = 2UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE simd_cdouble_t() noexcept : value( _mm_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE simd_cdouble_t( __m128d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t i ) const noexcept { return reinterpret_cast<const Type*>( &value )[i]; }
   __m128d value;  // Contains 1 64-bit double precision complex value
   enum : size_t { size = 1UL };
#else
   BLAZE_ALWAYS_INLINE simd_cdouble_t() noexcept : value( 0.0, 0.0 ) {}
   BLAZE_ALWAYS_INLINE simd_cdouble_t( Type v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE Type operator[]( size_t /*i*/ ) const noexcept { return value; }
   Type value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cdouble_t( const simd_cf64_t<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE simd_cdouble_t& operator=( const simd_cf64_t<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSIMDTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template<> struct IsSIMDType< simd_int8_t    > : public TrueType {};
template<> struct IsSIMDType< simd_uint8_t   > : public TrueType {};
template<> struct IsSIMDType< simd_cint8_t   > : public TrueType {};
template<> struct IsSIMDType< simd_cuint8_t  > : public TrueType {};

template<> struct IsSIMDType< simd_int16_t   > : public TrueType {};
template<> struct IsSIMDType< simd_uint16_t  > : public TrueType {};
template<> struct IsSIMDType< simd_cint16_t  > : public TrueType {};
template<> struct IsSIMDType< simd_cuint16_t > : public TrueType {};

template<> struct IsSIMDType< simd_int32_t   > : public TrueType {};
template<> struct IsSIMDType< simd_uint32_t  > : public TrueType {};
template<> struct IsSIMDType< simd_cint32_t  > : public TrueType {};
template<> struct IsSIMDType< simd_cuint32_t > : public TrueType {};

template<> struct IsSIMDType< simd_int64_t   > : public TrueType {};
template<> struct IsSIMDType< simd_uint64_t  > : public TrueType {};
template<> struct IsSIMDType< simd_cint64_t  > : public TrueType {};
template<> struct IsSIMDType< simd_cuint64_t > : public TrueType {};

template<> struct IsSIMDType< simd_float_t   > : public TrueType {};
template<> struct IsSIMDType< simd_cfloat_t  > : public TrueType {};

template<> struct IsSIMDType< simd_double_t  > : public TrueType {};
template<> struct IsSIMDType< simd_cdouble_t > : public TrueType {};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
