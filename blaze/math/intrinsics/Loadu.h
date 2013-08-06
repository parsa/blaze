//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Loadu.h
//  \brief Header file for the intrinsic unaligned load functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_LOADU_H_
#define _BLAZE_MATH_INTRINSICS_LOADU_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/Vectorization.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/Integral.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for intrinsic unaligned load operations.
// \ingroup intrinsics
//
// This helper structure provides the mapping between the size of an integral data type and the
// according intrinsic unaligned load function. Note that the type \a T must be an integral data
// type. Instantiating the Loadu class with a non-integral data type results in a compilation
// error.
*/
template< typename T  // Type of the integral
        , size_t N >  // Size of the integral
struct Loadu;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SPECIALIZATIONS OF THE LOADU CLASS TEMPLATE
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Loadu class template for 2-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Loadu<T,2UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int16_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type loadu( const T* address )
   {
#if BLAZE_AVX2_MODE
      return _mm256_loadu_si256( reinterpret_cast<const __m256i*>( address ) );
#elif BLAZE_SSE2_MODE
      return _mm_loadu_si128( reinterpret_cast<const __m128i*>( address ) );
#else
      return *address;
#endif
   }
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_INTEGRAL_TYPE( T );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Loadu class template for 4-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Loadu<T,4UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int32_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type loadu( const T* address )
   {
#if BLAZE_MIC_MODE
      __m512i v1;
      v1 = _mm512_loadunpacklo_epi32( v1, address );
      v1 = _mm512_loadunpackhi_epi32( v1, address+16UL );
      return v1;
#elif BLAZE_AVX2_MODE
      return _mm256_loadu_si256( reinterpret_cast<const __m256i*>( address ) );
#elif BLAZE_SSE2_MODE
      return _mm_loadu_si128( reinterpret_cast<const __m128i*>( address ) );
#else
      return *address;
#endif
   }
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_INTEGRAL_TYPE( T );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Loadu class template for 8-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Loadu<T,8UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int64_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type loadu( const T* address )
   {
#if BLAZE_MIC_MODE
      __m512i v1;
      v1 = _mm512_loadunpacklo_epi64( v1, address );
      v1 = _mm512_loadunpackhi_epi64( v1, address+8UL );
      return v1;
#elif BLAZE_AVX2_MODE
      return _mm256_loadu_si256( reinterpret_cast<const __m256i*>( address ) );
#elif BLAZE_SSE2_MODE
      return _mm_loadu_si128( reinterpret_cast<const __m128i*>( address ) );
#else
      return *address;
#endif
   }
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_INTEGRAL_TYPE( T );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INTRINSIC LOADU FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Loads a vector of integral values.
// \ingroup intrinsics
//
// \param address The first integral value to be loaded.
// \return The loaded vector of integral values.
//
// This function loads a vector of integral values. In contrast to the according load function,
// the given address is not required to be properly aligned.
*/
template< typename T >  // Type of the integral value
inline typename EnableIf< IsIntegral<T>, Loadu<T,sizeof(T)> >::Type::Type
   loadu( const T* address )
{
   return Loadu<T,sizeof(T)>::loadu( address );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'float' values.
// \ingroup intrinsics
//
// \param address The first 'float' value to be loaded.
// \return The loaded vector of 'float' values.
//
// This function loads a vector of 'float' values. In contrast to the according load function,
// the given address is not required to be properly aligned.
*/
inline sse_float_t loadu( const float* address )
{
#if BLAZE_MIC_MODE
   __m512 v1;
   v1 = _mm512_loadunpacklo_ps( v1, address );
   v1 = _mm512_loadunpackhi_ps( v1, address+16UL );
   return v1;
#elif BLAZE_AVX_MODE
   return _mm256_loadu_ps( address );
#elif BLAZE_SSE_MODE
   return _mm_loadu_ps( address );
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
//
// This function loads a vector of 'double' values. In contrast to the according load function,
// the given address is not required to be properly aligned.
*/
inline sse_double_t loadu( const double* address )
{
#if BLAZE_MIC_MODE
   __m512d v1;
   v1 = _mm512_loadunpacklo_pd( v1, address );
   v1 = _mm512_loadunpackhi_pd( v1, address+8UL );
   return v1;
#elif BLAZE_AVX_MODE
   return _mm256_loadu_pd( address );
#elif BLAZE_SSE2_MODE
   return _mm_loadu_pd( address );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'complex<float>' values.
// \ingroup intrinsics
//
// \param address The first 'complex<float>' value to be loaded.
// \return The loaded vector of 'complex<float>' values.
//
// This function loads a vector of 'complex<float>' values. In contrast to the according load
// function, the given address is not required to be properly aligned.
*/
inline sse_cfloat_t loadu( const complex<float>* address )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if BLAZE_AVX_MODE
   return _mm256_loadu_ps( reinterpret_cast<const float*>( address ) );
#elif BLAZE_SSE_MODE
   return _mm_loadu_ps( reinterpret_cast<const float*>( address ) );
#else
   return *address;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Loads a vector of 'complex<double>' values.
// \ingroup intrinsics
//
// \param address The first 'complex<double>' value to be loaded.
// \return The loaded vector of 'complex<double>' values.
//
// This function loads a vector of 'complex<double>' values. In contrast to the according load
// function, the given address is not required to be properly aligned.
*/
inline sse_cdouble_t loadu( const complex<double>* address )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if BLAZE_AVX_MODE
   return _mm256_loadu_pd( reinterpret_cast<const double*>( address ) );
#elif BLAZE_SSE2_MODE
   return _mm_loadu_pd( reinterpret_cast<const double*>( address ) );
#else
   return *address;
#endif
}
//*************************************************************************************************

} // namespace blaze

#endif
