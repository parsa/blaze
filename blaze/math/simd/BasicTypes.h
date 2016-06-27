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
/*\class blaze::SIMDPack
// \brief Base class for all SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD pack
struct SIMDPack {
   BLAZE_ALWAYS_INLINE T&       operator~() noexcept       { return *static_cast<T*>( this ); }
   BLAZE_ALWAYS_INLINE const T& operator~() const noexcept { return *static_cast<const T*>( this ); }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDi8
// \brief Base class for all 8-bit integral SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDi8 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDci8
// \brief Base class for all 8-bit integral complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDci8 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDi16
// \brief Base class for all 16-bit integral SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDi16 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDci16
// \brief Base class for all 16-bit integral complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDci16 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDi32
// \brief Base class for all 32-bit integral SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDi32 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDci32
// \brief Base class for all 32-bit integral complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDci32 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDi64
// \brief Base class for all 64-bit integral SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDi64 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDci64
// \brief Base class for all 64-bit integral complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDci64 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDf32
// \brief Base class for all single precision floating point SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDf32 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDcf32
// \brief Base class for all single precision floating point complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDcf32 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDf64
// \brief Base class for all double precision floating point SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDf64 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*\class blaze::SIMDcf64
// \brief Base class for all double precision floating point complex SIMD data types.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
template< typename T >  // Type of the SIMD element
struct SIMDcf64 : public SIMDPack< T >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  8-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDint8
// \brief SIMD type for 8-bit signed integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDint8 : public SIMDi8< SIMDint8 >
{
   using ValueType = int8_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDint8() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDint8( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDint8() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDint8( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 32 8-bit signed integral data values
   enum : size_t { size = 32UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDint8() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDint8( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 16 8-bit signed integral data values
   enum : size_t { size = 16UL };
#else
   BLAZE_ALWAYS_INLINE SIMDint8() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDint8( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDint8( const SIMDi8<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDint8& operator=( const SIMDi8<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::SIMDuint8
// \brief SIMD type for 8-bit unsigned integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDuint8 : public SIMDi8< SIMDuint8 >
{
   using ValueType = uint8_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDuint8() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDuint8( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDuint8() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint8( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 32 8-bit unsigned integral data values
   enum : size_t { size = 32UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDuint8() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint8( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 16 8-bit unsigned integral data values
   enum : size_t { size = 16UL };
#else
   BLAZE_ALWAYS_INLINE SIMDuint8() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDuint8( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDuint8( const SIMDi8<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDuint8& operator=( const SIMDi8<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  8-BIT INTEGRAL COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDcint8
// \brief SIMD type for 8-bit signed integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcint8 : public SIMDci8< SIMDcint8 >
{
   using ValueType = complex<int8_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcint8() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcint8( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDcint8() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint8( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 16 8-bit signed integral complex values
   enum : size_t { size = 16UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcint8() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint8( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 8 8-bit signed integral complex values
   enum : size_t { size = 8UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcint8() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcint8( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcint8( const SIMDci8<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcint8& operator=( const SIMDci8<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::SIMDcuint8
// \brief SIMD type for 8-bit unsigned integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcuint8 : public SIMDci8< SIMDcuint8 >
{
   using ValueType = complex<uint8_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint8() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint8( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint8() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint8( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 16 8-bit unsigned integral complex values
   enum : size_t { size = 16UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint8() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint8( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 8 8-bit unsigend integral complex values
   enum : size_t { size = 8UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcuint8() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint8( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcuint8( const SIMDci8<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcuint8& operator=( const SIMDci8<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  16-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDint16
// \brief SIMD type for 16-bit signed integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDint16 : public SIMDi16< SIMDint16 >
{
   using ValueType = int16_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDint16() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDint16( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDint16() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDint16( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 16 16-bit signed integral data values
   enum : size_t { size = 16UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDint16() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDint16( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 8 16-bit signed integral data values
   enum : size_t { size = 8UL };
#else
   BLAZE_ALWAYS_INLINE SIMDint16() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDint16( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDint16( const SIMDi16<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDint16& operator=( const SIMDi16<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::SIMDuint16
// \brief SIMD type for 16-bit unsigned integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDuint16 : public SIMDi16< SIMDuint16 >
{
   using ValueType = uint16_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDuint16() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDuint16( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDuint16() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint16( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 16 16-bit unsigned integral data values
   enum : size_t { size = 16UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDuint16() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint16( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 8 16-bit unsigned integral data values
   enum : size_t { size = 8UL };
#else
   BLAZE_ALWAYS_INLINE SIMDuint16() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDuint16( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDuint16( const SIMDi16<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDuint16& operator=( const SIMDi16<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  16-BIT INTEGRAL COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDcint16
// \brief SIMD type for 16-bit signed integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcint16 : public SIMDci16< SIMDcint16 >
{
   using ValueType = complex<int16_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcint16() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcint16( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDcint16() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint16( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 8 16-bit signed integral complex values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcint16() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint16( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 4 16-bit signed integral complex values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcint16() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcint16( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcint16( const SIMDci16<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcint16& operator=( const SIMDci16<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::SIMDcuint16
// \brief SIMD type for 16-bit unsigned integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcuint16 : public SIMDci16< SIMDcuint16 >
{
   using ValueType = complex<uint16_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint16() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint16( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint16() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint16( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 8 16-bit unsigned integral complex values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint16() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint16( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 4 16-bit unsigned integral complex values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcuint16() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint16( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcuint16( const SIMDci16<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcuint16& operator=( const SIMDci16<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDint32
// \brief SIMD type for 32-bit signed integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDint32 : public SIMDi32< SIMDint32 >
{
   using ValueType = int32_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDint32() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE SIMDint32( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512i value;  // Contains 16 32-bit signed integral data values
   enum : size_t { size = 16UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDint32() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDint32( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 8 32-bit signed integral data values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDint32() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDint32( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 4 32-bit signed integral data values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE SIMDint32() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDint32( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDint32( const SIMDi32<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDint32& operator=( const SIMDi32<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::SIMDuint32
// \brief SIMD type for 32-bit unsigned integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDuint32 : public SIMDi32< SIMDuint32 >
{
   using ValueType = uint32_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDuint32() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint32( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512i value;  // Contains 16 32-bit unsigned integral data values
   enum : size_t { size = 16UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDuint32() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint32( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 8 32-bit unsigned integral data values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDuint32() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint32( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 4 32-bit unsigned integral data values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE SIMDuint32() noexcept : value( 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDuint32( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDuint32( const SIMDi32<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDuint32& operator=( const SIMDi32<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  32-BIT INTEGRAL COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDcint32
// \brief SIMD type for 32-bit signed integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcint32 : public SIMDci32< SIMDcint32 >
{
   using ValueType = complex<int32_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcint32() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint32( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512i value;  // Contains 8 32-bit signed integral complex values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDcint32() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint32( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 4 32-bit signed integral complex values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcint32() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint32( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 2 32-bit signed integral complex values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcint32() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcint32( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcint32( const SIMDci32<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcint32& operator=( const SIMDci32<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::SIMDcuint32
// \brief SIMD type for 32-bit unsigned integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcuint32 : public SIMDci32< SIMDcuint32 >
{
   using ValueType = complex<uint32_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint32() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint32( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512i value;  // Contains 8 32-bit unsigned integral complex values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint32() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint32( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 4 32-bit unsigned integral complex values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint32() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint32( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 2 32-bit unsigned integral complex values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcuint32() noexcept : value( 0, 0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint32( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcuint32( const SIMDci32<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcuint32& operator=( const SIMDci32<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  64-BIT INTEGRAL SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDint64
// \brief SIMD type for 64-bit integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDint64 : public SIMDi64< SIMDint64 >
{
   using ValueType = int64_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDint64() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE SIMDint64( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512i value;  // Contains 8 64-bit signed integral data values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDint64() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDint64( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 4 64-bit signed integral data values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDint64() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDint64( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 2 64-bit signed integral data values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE SIMDint64() noexcept : value( 0L ) {}
   BLAZE_ALWAYS_INLINE SIMDint64( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDint64( const SIMDi64<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDint64& operator=( const SIMDi64<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::SIMDuint64
// \brief SIMD type for 64-bit unsigned integral data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDuint64 : public SIMDi64< SIMDuint64 >
{
   using ValueType = uint64_t;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDuint64() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint64( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512i value;  // Contains 8 64-bit unsigned integral data values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDuint64() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint64( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 4 64-bit unsigned integral data values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDuint64() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDuint64( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 2 64-bit unsigned integral data values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE SIMDuint64() noexcept : value( 0L ) {}
   BLAZE_ALWAYS_INLINE SIMDuint64( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDuint64( const SIMDi64<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDuint64& operator=( const SIMDi64<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  64-BIT INTEGRAL COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDcint64
// \brief SIMD type for 64-bit signed integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcint64 : public SIMDci64< SIMDcint64 >
{
   using ValueType = complex<int64_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcint64() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint64( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512i value;  // Contains 4 64-bit signed integral complex values
   enum : size_t { size = 4UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDcint64() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint64( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 2 64-bit signed integral complex values
   enum : size_t { size = 2UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcint64() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDcint64( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 1 64-bit signed integral complex values
   enum : size_t { size = 1UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcint64() noexcept : value( 0L, 0L ) {}
   BLAZE_ALWAYS_INLINE SIMDcint64( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcint64( const SIMDci64<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcint64& operator=( const SIMDci64<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::SIMDcuint64
// \brief SIMD type for 64-bit unsigned integral complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcuint64 : public SIMDci64< SIMDcuint64 >
{
   using ValueType = complex<uint64_t>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint64() noexcept : value( _mm512_setzero_epi32() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint64( __m512i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512i value;  // Contains 4 64-bit unsigned integral complex values
   enum : size_t { size = 4UL };
#elif BLAZE_AVX2_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint64() noexcept : value( _mm256_setzero_si256() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint64( __m256i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256i value;  // Contains 2 64-bit unsigned integral complex values
   enum : size_t { size = 2UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcuint64() noexcept : value( _mm_setzero_si128() ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint64( __m128i v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128i value;  // Contains 1 64-bit unsigned integral complex values
   enum : size_t { size = 1UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcuint64() noexcept : value( 0L, 0L ) {}
   BLAZE_ALWAYS_INLINE SIMDcuint64( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcuint64( const SIMDci64<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcuint64& operator=( const SIMDci64<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SINGLE PRECISION FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDfloat
// \brief SIMD type for 32-bit single precision floating point data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDfloat : public SIMDf32< SIMDfloat >
{
   using ValueType = float;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDfloat() noexcept : value( _mm512_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE SIMDfloat( __m512 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512 value;  // Contains 16 32-bit single precision floating point values
   enum : size_t { size = 16UL };
#elif BLAZE_AVX_MODE
   BLAZE_ALWAYS_INLINE SIMDfloat() noexcept : value( _mm256_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE SIMDfloat( __m256 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256 value;  // Contains 8 32-bit single precision floating point values
   enum : size_t { size = 8UL };
#elif BLAZE_SSE_MODE
   BLAZE_ALWAYS_INLINE SIMDfloat() noexcept : value( _mm_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE SIMDfloat( __m128 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128 value;  // Contains 4 32-bit single precision floating point values
   enum : size_t { size = 4UL };
#else
   BLAZE_ALWAYS_INLINE SIMDfloat() noexcept : value( 0.0F ) {}
   BLAZE_ALWAYS_INLINE SIMDfloat( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDfloat( const SIMDf32<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDfloat& operator=( const SIMDf32<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SINGLE PRECISION FLOATING POINT COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDcfloat
// \brief SIMD type for 32-bit single precision complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcfloat : public SIMDcf32< SIMDcfloat >
{
   using ValueType = complex<float>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcfloat() noexcept : value( _mm512_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE SIMDcfloat( __m512 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512 value;  // Contains 8 32-bit single precision complex values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX_MODE
   BLAZE_ALWAYS_INLINE SIMDcfloat() noexcept : value( _mm256_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE SIMDcfloat( __m256 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256 value;  // Contains 4 32-bit single precision complex values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE_MODE
   BLAZE_ALWAYS_INLINE SIMDcfloat() noexcept : value( _mm_setzero_ps() ) {}
   BLAZE_ALWAYS_INLINE SIMDcfloat( __m128 v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128 value;  // Contains 2 32-bit single precision complex values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcfloat() noexcept : value( 0.0F, 0.0F ) {}
   BLAZE_ALWAYS_INLINE SIMDcfloat( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcfloat( const SIMDcf32<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcfloat& operator=( const SIMDcf32<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DOUBLE PRECISION FLOATING POINT SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDdouble
// \brief SIMD type for 64-bit double precision floating point data values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDdouble : public SIMDf64< SIMDdouble >
{
   using ValueType = double;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDdouble() noexcept : value( _mm512_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE SIMDdouble( __m512d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512d value;  // Contains 8 64-bit double precision floating point values
   enum : size_t { size = 8UL };
#elif BLAZE_AVX_MODE
   BLAZE_ALWAYS_INLINE SIMDdouble() noexcept : value( _mm256_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE SIMDdouble( __m256d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256d value;  // Contains 4 64-bit double precision floating point values
   enum : size_t { size = 4UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDdouble() noexcept : value( _mm_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE SIMDdouble( __m128d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128d value;  // Contains 2 64-bit double precision floating point values
   enum : size_t { size = 2UL };
#else
   BLAZE_ALWAYS_INLINE SIMDdouble() noexcept : value( 0.0 ) {}
   BLAZE_ALWAYS_INLINE SIMDdouble( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDdouble( const SIMDf64<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDdouble& operator=( const SIMDf64<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DOUBLE PRECISION FLOATING POINT COMPLEX SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::SIMDcdouble
// \brief SIMD type for 64-bit double precision complex values.
// \ingroup simd
*/
/*! \cond BLAZE_INTERNAL */
struct SIMDcdouble : public SIMDcf64< SIMDcdouble >
{
   using ValueType = complex<double>;

#if BLAZE_MIC_MODE
   BLAZE_ALWAYS_INLINE SIMDcdouble() noexcept : value( _mm512_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE SIMDcdouble( __m512d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m512d value;  // Contains 4 64-bit double precision complex value
   enum : size_t { size = 4UL };
#elif BLAZE_AVX_MODE
   BLAZE_ALWAYS_INLINE SIMDcdouble() noexcept : value( _mm256_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE SIMDcdouble( __m256d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m256d value;  // Contains 2 64-bit double precision complex value
   enum : size_t { size = 2UL };
#elif BLAZE_SSE2_MODE
   BLAZE_ALWAYS_INLINE SIMDcdouble() noexcept : value( _mm_setzero_pd() ) {}
   BLAZE_ALWAYS_INLINE SIMDcdouble( __m128d v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t i ) const noexcept { return reinterpret_cast<const ValueType*>( &value )[i]; }
   __m128d value;  // Contains 1 64-bit double precision complex value
   enum : size_t { size = 1UL };
#else
   BLAZE_ALWAYS_INLINE SIMDcdouble() noexcept : value( 0.0, 0.0 ) {}
   BLAZE_ALWAYS_INLINE SIMDcdouble( ValueType v ) noexcept : value( v ) {}
   BLAZE_ALWAYS_INLINE ValueType operator[]( size_t /*i*/ ) const noexcept { return value; }
   ValueType value;
   enum : size_t { size = 1UL };
#endif

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcdouble( const SIMDcf64<T>& v ) noexcept : value( (~v).value ) {}

   template< typename T >
   BLAZE_ALWAYS_INLINE SIMDcdouble& operator=( const SIMDcf64<T>& v ) noexcept { value = (~v).value; return *this; }
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIMD OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of two SIMD packs.
// \ingroup simd
//
// \param lhs The left-hand side SIMD operand for the addition.
// \param rhs The right-hand side SIMD operand for the addition.
// \return Reference to the left-hand side SIMD operand.
*/
template< typename T1    // Type of the left-hand side SIMD operand
        , typename T2 >  // Type of the right-hand side SIMD operand
BLAZE_ALWAYS_INLINE T1& operator+=( SIMDPack<T1>& lhs, const SIMDPack<T2>& rhs )
{
   (~lhs) = (~lhs) + (~rhs);
   return ~lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of two SIMD packs.
// \ingroup simd
//
// \param lhs The left-hand side SIMD operand for the subtraction.
// \param rhs The right-hand side SIMD operand for the subtraction.
// \return Reference to the left-hand side SIMD operand.
*/
template< typename T1    // Type of the left-hand side SIMD operand
        , typename T2 >  // Type of the right-hand side SIMD operand
BLAZE_ALWAYS_INLINE T1& operator-=( SIMDPack<T1>& lhs, const SIMDPack<T2>& rhs )
{
   (~lhs) = (~lhs) - (~rhs);
   return ~lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of two SIMD packs.
// \ingroup simd
//
// \param lhs The left-hand side SIMD operand for the multiplication.
// \param rhs The right-hand side SIMD operand for the multiplication.
// \return Reference to the left-hand side SIMD operand.
*/
template< typename T1    // Type of the left-hand side SIMD operand
        , typename T2 >  // Type of the right-hand side SIMD operand
BLAZE_ALWAYS_INLINE T1& operator*=( SIMDPack<T1>& lhs, const SIMDPack<T2>& rhs )
{
   (~lhs) = (~lhs) * (~rhs);
   return ~lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of two SIMD packs.
// \ingroup simd
//
// \param lhs The left-hand side SIMD operand for the division.
// \param rhs The right-hand side SIMD operand for the division.
// \return Reference to the left-hand side SIMD operand.
*/
template< typename T1    // Type of the left-hand side SIMD operand
        , typename T2 >  // Type of the right-hand side SIMD operand
BLAZE_ALWAYS_INLINE T1& operator/=( SIMDPack<T1>& lhs, const SIMDPack<T2>& rhs )
{
   (~lhs) = (~lhs) / (~rhs);
   return ~lhs;
}
//*************************************************************************************************




//=================================================================================================
//
//  ISSIMDTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template<> struct IsSIMDType< SIMDint8    > : public TrueType {};
template<> struct IsSIMDType< SIMDuint8   > : public TrueType {};
template<> struct IsSIMDType< SIMDcint8   > : public TrueType {};
template<> struct IsSIMDType< SIMDcuint8  > : public TrueType {};

template<> struct IsSIMDType< SIMDint16   > : public TrueType {};
template<> struct IsSIMDType< SIMDuint16  > : public TrueType {};
template<> struct IsSIMDType< SIMDcint16  > : public TrueType {};
template<> struct IsSIMDType< SIMDcuint16 > : public TrueType {};

template<> struct IsSIMDType< SIMDint32   > : public TrueType {};
template<> struct IsSIMDType< SIMDuint32  > : public TrueType {};
template<> struct IsSIMDType< SIMDcint32  > : public TrueType {};
template<> struct IsSIMDType< SIMDcuint32 > : public TrueType {};

template<> struct IsSIMDType< SIMDint64   > : public TrueType {};
template<> struct IsSIMDType< SIMDuint64  > : public TrueType {};
template<> struct IsSIMDType< SIMDcint64  > : public TrueType {};
template<> struct IsSIMDType< SIMDcuint64 > : public TrueType {};

template<> struct IsSIMDType< SIMDfloat   > : public TrueType {};
template<> struct IsSIMDType< SIMDcfloat  > : public TrueType {};

template<> struct IsSIMDType< SIMDdouble  > : public TrueType {};
template<> struct IsSIMDType< SIMDcdouble > : public TrueType {};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
