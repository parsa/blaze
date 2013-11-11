//=================================================================================================
/*!
//  \file blaze/math/intrinsics/Set.h
//  \brief Header file for the intrinsic set functionality
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
#include <blaze/util/StaticAssert.h>


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
#if BLAZE_AVX2_MODE
      return _mm256_set1_epi16( value );
#elif BLAZE_SSE2_MODE
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
#elif BLAZE_AVX2_MODE
      return _mm256_set1_epi32( value );
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
#elif BLAZE_AVX2_MODE
      return _mm256_set1_epi64x( value );
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


//*************************************************************************************************
/*!\brief Sets all values in the vector to the given 'complex<float>' value.
// \ingroup intrinsics
//
// \param value The given 'complex<float>' value.
// \return The set vector of 'complex<float>' values.
*/
inline sse_cfloat_t set( const complex<float>& value )
{
#if BLAZE_MIC_MODE
   return _mm512_set_ps( value.imag(), value.real(), value.imag(), value.real(),
                         value.imag(), value.real(), value.imag(), value.real(),
                         value.imag(), value.real(), value.imag(), value.real(),
                         value.imag(), value.real(), value.imag(), value.real() );
#elif BLAZE_AVX_MODE
   return _mm256_set_ps( value.imag(), value.real(), value.imag(), value.real(),
                         value.imag(), value.real(), value.imag(), value.real() );
#elif BLAZE_SSE_MODE
   return _mm_set_ps( value.imag(), value.real(), value.imag(), value.real() );
#else
   return value;
#endif
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets all values in the vector to the given 'complex<double>' value.
// \ingroup intrinsics
//
// \param value The given 'complex<double>' value.
// \return The set vector of 'complex<double>' values.
*/
inline sse_cdouble_t set( const complex<double>& value )
{
#if BLAZE_MIC_MODE
   return _mm512_set_pd( value.imag(), value.real(), value.imag(), value.real(),
                         value.imag(), value.real(), value.imag(), value.real() );
#elif BLAZE_AVX_MODE
   return _mm256_set_pd( value.imag(), value.real(), value.imag(), value.real() );
#elif BLAZE_SSE2_MODE
   return _mm_set_pd( value.imag(), value.real() );
#else
   return value;
#endif
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
