//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Storeu.h
//  \brief Header file for the intrinsic unaligned store functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_STOREU_H_
#define _BLAZE_MATH_INTRINSICS_STOREU_H_


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
/*!\brief Auxiliary helper struct for intrinsic aligned store operations.
// \ingroup intrinsics
//
// This helper structure provides the mapping between the size of an integral data type and the
// according intrinsic aligned store function. Note that the type \a T must be an integral data
// type. Instantiating the Storeu class with a non-integral data type results in a compilation
// error.
*/
template< typename T  // Type of the integral
        , size_t N >  // Size of the integral
struct Storeu;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SPECIALIZATIONS OF THE STOREU CLASS TEMPLATE
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Storeu class template for 2-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Storeu<T,2UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int16_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline void storeu( T* address, const Type& value )
   {
#if BLAZE_AVX2_MODE
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
      _mm_storeu_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
      *address = value.value;
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
/*!\brief Specialization of the Storeu class template for 4-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Storeu<T,4UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int32_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline void storeu( T* address, const Type& value )
   {
#if BLAZE_MIC_MODE
      _mm512_packstorelo_epi32( address, value.value );
      _mm512_packstorehi_epi32( address+16UL, value.value );
#elif BLAZE_AVX2_MODE
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
      _mm_storeu_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
      *address = value.value;
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
/*!\brief Specialization of the Storeu class template for 8-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Storeu<T,8UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int64_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline void storeu( T* address, const Type& value )
   {
#if BLAZE_MIC_MODE
      _mm512_packstorelo_epi64( address, value.value );
      _mm512_packstorehi_epi64( address+8UL, value.value );
#elif BLAZE_AVX2_MODE
      _mm256_storeu_si256( reinterpret_cast<__m256i*>( address ), value.value );
#elif BLAZE_SSE2_MODE
      _mm_storeu_si128( reinterpret_cast<__m128i*>( address ), value.value );
#else
      *address = value.value;
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
//  INTRINSIC STOREU FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Unaligned store of a vector of integral values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The integral vector to be stored.
// \return void
//
// This function stores a vector of integral values. In contrast to the according store function,
// the given address is not required to be properly aligned.
*/
template< typename T >  // Type of the integral value
inline typename EnableIf< IsIntegral<T> >::Type
   storeu( T* address, const typename Storeu<T,sizeof(T)>::Type& value )
{
   Storeu<T,sizeof(T)>::storeu( address, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of a vector of 'float' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'float' vector to be stored.
// \return void
//
// This function stores a vector of 'float' values. In contrast to the according store function,
// the given address is not required to be properly aligned.
*/
inline void storeu( float* address, const sse_float_t& value )
{
#if BLAZE_MIC_MODE
   _mm512_packstorelo_ps( address     , value.value );
   _mm512_packstorehi_ps( address+16UL, value.value );
#elif BLAZE_AVX_MODE
   _mm256_storeu_ps( address, value.value );
#elif BLAZE_SSE_MODE
   _mm_storeu_ps( address, value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of a vector of 'double' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'double' vector to be stored.
// \return void
//
// This function stores a vector of 'double' values. In contrast to the according store function,
// the given address is not required to be properly aligned.
*/
inline void storeu( double* address, const sse_double_t& value )
{
#if BLAZE_MIC_MODE
   _mm512_packstorelo_pd( address    , value.value );
   _mm512_packstorehi_pd( address+8UL, value.value );
#elif BLAZE_AVX_MODE
   _mm256_storeu_pd( address, value.value );
#elif BLAZE_SSE2_MODE
   _mm_storeu_pd( address, value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of a vector of 'complex<float>' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'complex<float>' vector to be stored.
// \return void
//
// This function stores a vector of 'complex<float>' values. In contrast to the according store
// function, the given address is not required to be properly aligned.
*/
inline void storeu( complex<float>* address, const sse_cfloat_t& value )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if BLAZE_MIC_MODE
   _mm512_packstorelo_ps( reinterpret_cast<float*>( address     ), value.value );
   _mm512_packstorehi_ps( reinterpret_cast<float*>( address+8UL ), value.value );
#elif BLAZE_AVX_MODE
   _mm256_storeu_ps( reinterpret_cast<float*>( address ), value.value );
#elif BLAZE_SSE_MODE
   _mm_storeu_ps( reinterpret_cast<float*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of a vector of 'complex<double>' values.
// \ingroup intrinsics
//
// \param address The target address.
// \param value The 'complex<double>' vector to be stored.
// \return void
//
// This function stores a vector of 'complex<double>' values. In contrast to the according store
// function, the given address is not required to be properly aligned.
*/
inline void storeu( complex<double>* address, const sse_cdouble_t& value )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if BLAZE_MIC_MODE
   _mm512_packstorelo_pd( reinterpret_cast<double*>( address     ), value.value );
   _mm512_packstorehi_pd( reinterpret_cast<double*>( address+4UL ), value.value );
#elif BLAZE_AVX_MODE
   _mm256_storeu_pd( reinterpret_cast<double*>( address ), value.value );
#elif BLAZE_SSE2_MODE
   _mm_storeu_pd( reinterpret_cast<double*>( address ), value.value );
#else
   *address = value.value;
#endif
}
//*************************************************************************************************

} // namespace blaze

#endif
