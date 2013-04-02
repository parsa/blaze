//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Set.h
//  \brief Header file for the intrinsic set functionality
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

#ifndef _BLAZE_MATH_INTRINSICS_SET_H_
#define _BLAZE_MATH_INTRINSICS_SET_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/intrinsics/BasicTypes.h>
#include <blaze/system/Vectorization.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Integral.h>
#include <blaze/util/EnableIf.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for intrinsic set operations.
// \ingroup intrinsics
//
// This helper structure provides the mapping between the size of an integral data type and the
// according intrinsic set function. Note that the type \a T must be an integral data type.
// Instantiating the Set class with a non-integral data type results in a compilation error.
*/
template< typename T  // Type of the integral
        , size_t N >  // Size of the integral
struct Set;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SPECIALIZATIONS OF THE SET CLASS TEMPLATE
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Set class template for 2-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Set<T,2UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int16_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type set( T value )
   {
#if BLAZE_SSE2_MODE
      return _mm_set1_epi16( value );
#else
      return value;
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
/*!\brief Specialization of the Set class template for 4-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Set<T,4UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int32_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type set( T value )
   {
#if BLAZE_MIC_MODE
      return _mm512_set1_epi32( value );
#elif BLAZE_SSE2_MODE
      return _mm_set1_epi32( value );
#else
      return value;
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
/*!\brief Specialization of the Set class template for 8-byte integral data types.
// \ingroup intrinsics
*/
template< typename T >  // Type of the integral
struct Set<T,8UL>
{
 public:
   //**Type definitions****************************************************************************
   typedef sse_int64_t  Type;
   //**********************************************************************************************

   //**Set function********************************************************************************
   static inline Type set( T value )
   {
#if BLAZE_MIC_MODE
      return _mm512_set1_epi64( value );
#elif BLAZE_SSE2_MODE
      return _mm_set1_epi64( value );
#else
      return value;
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
//  INTRINSIC SET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Sets all values in the vector to the given integral value.
// \ingroup intrinsics
//
// \param value The given integral value.
// \return The set vector of integral values.
*/
template< typename T >  // Type of the integral value
inline typename EnableIf< IsIntegral<T>, Set<T,sizeof(T)> >::Type::Type
   set( T value )
{
   return Set<T,sizeof(T)>::set( value );
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

} // namespace blaze

#endif
