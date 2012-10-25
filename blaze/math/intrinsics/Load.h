//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Load.h
//  \brief Header file for the intrinsic load functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_LOAD_H_
#define _BLAZE_MATH_INTRINSICS_LOAD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/SSE.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Integral.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for intrinsic load operations.
// \ingroup intrinsics
//
// This helper structure provides the mapping between the size of an integral data type and the
// according intrinsic load function. Note that the type \a T must be an integral data type.
// Instantiating the Load class with a non-integral data type results in a compilation error.
*/
template< typename T  // Type of the integral
        , size_t N >  // Size of the integral
struct Load;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SPECIALIZATIONS OF THE LOAD CLASS TEMPLATE
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Load class template for 2-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Load<T,2UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int16_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type load( const T* address )
   {
#if BLAZE_SSE2_MODE
      BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
      return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
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
/*!\brief Specialization of the Load class template for 4-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Load<T,4UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int32_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type load( const T* address )
   {
#if BLAZE_MIC_MODE
      BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 64UL ), "Invalid alignment detected" );
      return _mm512_load_epi32( address );
#elif BLAZE_SSE2_MODE
      BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
      return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
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
/*!\brief Specialization of the Load class template for 8-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Load<T,8UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int64_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type load( const T* address )
   {
#if BLAZE_MIC_MODE
      BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 64UL ), "Invalid alignment detected" );
      return _mm512_load_epi64( address );
#elif BLAZE_SSE2_MODE
      BLAZE_INTERNAL_ASSERT( !( reinterpret_cast<size_t>( address ) % 16UL ), "Invalid alignment detected" );
      return _mm_load_si128( reinterpret_cast<const __m128i*>( address ) );
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
//  INTRINSIC LOAD FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Loads a vector of integral values.
// \ingroup intrinsics
//
// \param address The first integral value to be loaded.
// \return The loaded vector of integral values.
*/
template< typename T >  // Type of the integral value
inline typename EnableIf< IsIntegral<T>, Load<T,sizeof(T)> >::Type::Type
   load( const T* address )
{
   return Load<T,sizeof(T)>::load( address );
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

} // namespace blaze

#endif
