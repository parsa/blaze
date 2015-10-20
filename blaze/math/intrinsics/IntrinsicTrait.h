//=================================================================================================
/*!
//  \file blaze/math/intrinsics/IntrinsicTrait.h
//  \brief Header file for the intrinsic trait
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

#ifndef _BLAZE_MATH_INTRINSICS_INTRINSICTRAIT_H_
#define _BLAZE_MATH_INTRINSICS_INTRINSICTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/Vectorization.h>
#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blaze/util/typetraits/RemoveCV.h>


namespace blaze {

//=================================================================================================
//
//  CLASS INTRINSICTRAITHELPER
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IntrinsicTrait class template.
// \ingroup intrinsics
//
// This helper structure provides the mapping between the size of an integral data type and the
// according intrinsic data type.
*/
template< bool C, size_t N >
struct IntrinsicTraitHelper;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief IntrinsicTraitHelper specialization for 1-byte integral data types.
// \ingroup intrinsics
*/
#if BLAZE_AVX2_MODE
template<>
struct IntrinsicTraitHelper<false,1UL>
{
   typedef simd_int8_t  Type;
   enum { size           = 32,
          addition       = 1,
          subtraction    = 1,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 1,
          conjugate      = 1 };
};
#else
template<>
struct IntrinsicTraitHelper<false,1UL>
{
   typedef simd_int8_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 16 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = 0,
          division       = 0,
          absoluteValue  = BLAZE_SSSE3_MODE,
          conjugate      = 1 };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief IntrinsicTraitHelper specialization for 2-byte integral data types.
// \ingroup intrinsics
*/
#if BLAZE_AVX2_MODE
template<>
struct IntrinsicTraitHelper<false,2UL>
{
   typedef simd_int16_t  Type;
   enum { size           = 16,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 0,
          absoluteValue  = 1,
          conjugate      = 1 };
};
#else
template<>
struct IntrinsicTraitHelper<false,2UL>
{
   typedef simd_int16_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 8 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = BLAZE_SSE2_MODE,
          division       = 0,
          absoluteValue  = BLAZE_SSSE3_MODE,
          conjugate      = 1 };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief IntrinsicTraitHelper specialization for 4-byte integral data types.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitHelper<false,4UL>
{
   typedef simd_int32_t  Type;
   enum { size           = 16,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 1,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#elif BLAZE_AVX2_MODE
template<>
struct IntrinsicTraitHelper<false,4UL>
{
   typedef simd_int32_t  Type;
   enum { size           = 8,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 0,
          absoluteValue  = 1,
          conjugate      = 1 };
};
#else
template<>
struct IntrinsicTraitHelper<false,4UL>
{
   typedef simd_int32_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 4 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = BLAZE_SSE4_MODE,
          division       = 0,
          absoluteValue  = BLAZE_SSSE3_MODE,
          conjugate      = 1 };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief IntrinsicTraitHelper specialization for 8-byte integral data types.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitHelper<false,8UL>
{
   typedef simd_int64_t  Type;
   enum { size           = 8,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 1,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#elif BLAZE_AVX2_MODE
template<>
struct IntrinsicTraitHelper<false,8UL>
{
   typedef simd_int64_t  Type;
   enum { size           = 4,
          addition       = 1,
          subtraction    = 1,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#else
template<>
struct IntrinsicTraitHelper<false,8UL>
{
   typedef simd_int64_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 2 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief IntrinsicTraitHelper specialization for 1-byte integral complex data types.
// \ingroup intrinsics
*/
#if BLAZE_AVX2_MODE
template<>
struct IntrinsicTraitHelper<true,1UL>
{
   typedef simd_cint8_t  Type;
   enum { size           = 16,
          addition       = 1,
          subtraction    = 1,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 0 };
};
#else
template<>
struct IntrinsicTraitHelper<true,1UL>
{
   typedef simd_cint8_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 8 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 0 };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief IntrinsicTraitHelper specialization for 2-byte integral complex data types.
// \ingroup intrinsics
*/
#if BLAZE_AVX2_MODE
template<>
struct IntrinsicTraitHelper<true,2UL>
{
   typedef simd_cint16_t  Type;
   enum { size           = 8,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#else
template<>
struct IntrinsicTraitHelper<true,2UL>
{
   typedef simd_cint16_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 4 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = BLAZE_SSE2_MODE,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = BLAZE_SSE2_MODE };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief IntrinsicTraitHelper specialization for 4-byte integral complex data types.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitHelper<true,4UL>
{
   typedef simd_cint32_t  Type;
   enum { size           = 8,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#elif BLAZE_AVX2_MODE
template<>
struct IntrinsicTraitHelper<true,4UL>
{
   typedef simd_cint32_t  Type;
   enum { size           = 4,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#else
template<>
struct IntrinsicTraitHelper<true,4UL>
{
   typedef simd_cint32_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 2 ):( 1 ),
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = BLAZE_SSE4_MODE,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = BLAZE_SSE4_MODE };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief IntrinsicTraitHelper specialization for 8-byte integral complex data types.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitHelper<true,8UL>
{
   typedef simd_cint64_t  Type;
   enum { size           = 4,
          addition       = 1,
          subtraction    = 1,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 0 };
};
#elif BLAZE_AVX2_MODE
template<>
struct IntrinsicTraitHelper<true,8UL>
{
   typedef simd_cint64_t  Type;
   enum { size           = 2,
          addition       = 1,
          subtraction    = 1,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 0 };
};
#else
template<>
struct IntrinsicTraitHelper<true,8UL>
{
   typedef simd_cint64_t  Type;
   enum { size           = 1,
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 0 };
};
#endif
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS INTRINSICTRAITBASE
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Base template for the IntrinsicTraitBase class.
// \ingroup intrinsics
*/
template< typename T >
struct IntrinsicTraitBase
{
   typedef T  Type;
   enum { size           = 1,
          alignment      = AlignmentOf<T>::value,
          addition       = 0,
          subtraction    = 0,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 0 };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'char'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<char>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(char)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<char>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'signed char'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<signed char>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(signed char)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<signed char>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'unsigned char'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<unsigned char>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(unsigned char)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<unsigned char>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = 0,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'wchar_t'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<wchar_t>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(wchar_t)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<wchar_t>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'short'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<short>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(short)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<short>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'unsigned short'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<unsigned short>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(unsigned short)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<unsigned short>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = 0,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'int'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<int>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(int)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<int>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'unsigned int'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<unsigned int>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(unsigned int)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<unsigned int>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = 0,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'long'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<long>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(long)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<long>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'unsigned long'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase<unsigned long>
{
 private:
   typedef IntrinsicTraitHelper<false,sizeof(unsigned long)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf<unsigned long>::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = 0,
          conjugate      = Helper::conjugate };
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'float'.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
struct IntrinsicTraitBase<float>
{
   typedef simd_float_t  Type;
   enum { size           = ( 64UL / sizeof(float) ),
          alignment      = AlignmentOf<float>::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 1,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#elif BLAZE_AVX_MODE
template<>
struct IntrinsicTraitBase<float>
{
   typedef simd_float_t  Type;
   enum { size           = ( 32UL / sizeof(float) ),
          alignment      = AlignmentOf<float>::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 1,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#else
template<>
struct IntrinsicTraitBase<float>
{
   typedef simd_float_t  Type;
   enum { size           = ( BLAZE_SSE_MODE )?( 16UL / sizeof(float) ):( 1 ),
          alignment      = AlignmentOf<float>::value,
          addition       = BLAZE_SSE_MODE,
          subtraction    = BLAZE_SSE_MODE,
          multiplication = BLAZE_SSE_MODE,
          division       = BLAZE_SSE_MODE,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'double'.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitBase<double>
{
   typedef simd_double_t  Type;
   enum { size           = ( 64UL / sizeof(double) ),
          alignment      = AlignmentOf<double>::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 1,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#elif BLAZE_AVX_MODE
template<>
struct IntrinsicTraitBase<double>
{
   typedef simd_double_t  Type;
   enum { size           = ( 32UL / sizeof(double) ),
          alignment      = AlignmentOf<double>::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 1,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#else
template<>
struct IntrinsicTraitBase<double>
{
   typedef simd_double_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 16UL / sizeof(double) ):( 1 ),
          alignment      = AlignmentOf<double>::value,
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = BLAZE_SSE2_MODE,
          division       = BLAZE_SSE2_MODE,
          absoluteValue  = 0,
          conjugate      = 1 };
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<char>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<char> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(char)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<char> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<char> ) == 2UL*sizeof( char ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<signed char>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<signed char> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(signed char)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<signed char> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<signed char> ) == 2UL*sizeof( signed char ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<unsigned char>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<unsigned char> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(unsigned char)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<unsigned char> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = 0,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<unsigned char> ) == 2UL*sizeof( unsigned char ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<wchar_t>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<wchar_t> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(wchar_t)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<wchar_t> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<wchar_t> ) == 2UL*sizeof( wchar_t ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<short>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<short> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(short)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<short> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<short> ) == 2UL*sizeof( short ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<unsigned short>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<unsigned short> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(unsigned short)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<unsigned short> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = 0,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<unsigned short> ) == 2UL*sizeof( unsigned short ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<int>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<int> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(int)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<int> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<int> ) == 2UL*sizeof( int ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<unsigned int>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<unsigned int> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(unsigned int)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<unsigned int> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = 0,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<unsigned int> ) == 2UL*sizeof( unsigned int ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<long>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<long> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(long)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<long> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = Helper::absoluteValue,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<long> ) == 2UL*sizeof( long ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<unsigned long>'.
// \ingroup intrinsics
*/
template<>
struct IntrinsicTraitBase< complex<unsigned long> >
{
 private:
   typedef IntrinsicTraitHelper<true,sizeof(unsigned long)>  Helper;

 public:
   typedef Helper::Type  Type;
   enum { size           = Helper::size,
          alignment      = AlignmentOf< complex<unsigned long> >::value,
          addition       = Helper::addition,
          subtraction    = Helper::subtraction,
          multiplication = Helper::multiplication,
          division       = Helper::division,
          absoluteValue  = 0,
          conjugate      = Helper::conjugate };

   BLAZE_STATIC_ASSERT( sizeof( complex<unsigned long> ) == 2UL*sizeof( unsigned long ) );
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<float>'.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitBase< complex<float> >
{
   typedef simd_cfloat_t  Type;
   enum { size           = ( 64UL / sizeof(complex<float>) ),
          alignment      = AlignmentOf< complex<float> >::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };

   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );
};
#elif BLAZE_AVX_MODE
template<>
struct IntrinsicTraitBase< complex<float> >
{
   typedef simd_cfloat_t  Type;
   enum { size           = ( 32UL / sizeof(complex<float>) ),
          alignment      = AlignmentOf< complex<float> >::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };

   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );
};
#else
template<>
struct IntrinsicTraitBase< complex<float> >
{
   typedef simd_cfloat_t  Type;
   enum { size           = ( BLAZE_SSE_MODE )?( 16UL / sizeof(complex<float>) ):( 1 ),
          alignment      = AlignmentOf< complex<float> >::value,
          addition       = BLAZE_SSE_MODE,
          subtraction    = BLAZE_SSE_MODE,
          multiplication = BLAZE_SSE3_MODE,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = BLAZE_SSE_MODE };

   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );
};
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IntrinsicTraitBase class template for 'complex<double>'.
// \ingroup intrinsics
*/
#if BLAZE_MIC_MODE
template<>
struct IntrinsicTraitBase< complex<double> >
{
   typedef simd_cdouble_t  Type;
   enum { size           = ( 64UL / sizeof(complex<double>) ),
          alignment      = AlignmentOf< complex<double> >::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 0,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };

   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );
};
#elif BLAZE_AVX_MODE
template<>
struct IntrinsicTraitBase< complex<double> >
{
   typedef simd_cdouble_t  Type;
   enum { size           = ( 32UL / sizeof(complex<double>) ),
          alignment      = AlignmentOf< complex<double> >::value,
          addition       = 1,
          subtraction    = 1,
          multiplication = 1,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = 1 };

   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );
};
#else
template<>
struct IntrinsicTraitBase< complex<double> >
{
   typedef simd_cdouble_t  Type;
   enum { size           = ( BLAZE_SSE2_MODE )?( 16UL / sizeof(complex<double>) ):( 1 ),
          alignment      = AlignmentOf< complex<double> >::value,
          addition       = BLAZE_SSE2_MODE,
          subtraction    = BLAZE_SSE2_MODE,
          multiplication = BLAZE_SSE3_MODE,
          division       = 0,
          absoluteValue  = 0,
          conjugate      = BLAZE_SSE2_MODE };

   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );
};
#endif
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS INTRINSICTRAIT
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Intrinsic characteristics of data types.
// \ingroup intrinsics
//
// The IntrinsicTrait class template provides the intrinsic characteristics of a specific data
// type:
//
//  - The nested data type \a Type corresponds to the according packed, intrinsic data type. In
//    case the data type doesn't have an intrinsic representation, \a Type corresonds to the given
//    data type itself.
//  - The \a size value corresponds to the number of values of the given data type that are packed
//    together in one intrinsic vector type. In case the data type cannot be vectorized, \a size
//    is set to 1.
//  - If the data type can be involved in vectorized additions, the \a addition value is set to 1.
//    Otherwise, \a addition is set to 0.
//  - In case the data type supports vectorized subtractions, the \a subtraction value is set to 1.
//    Else it is set to 0.
//  - If the data type supports vectorized multiplications, the \a multiplication value is set to
//    1. If it cannot be used in multiplications, it is set to 0.
*/
template< typename T >
class IntrinsicTrait : public IntrinsicTraitBase< typename RemoveCV<T>::Type >
{};
//*************************************************************************************************

} // namespace blaze

#endif
