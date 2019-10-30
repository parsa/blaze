//=================================================================================================
/*!
//  \file blaze/math/functors/Clamp.h
//  \brief Header file for the Clamp functor
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_FUNCTORS_CLAMP_H_
#define _BLAZE_MATH_FUNCTORS_CLAMP_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/SIMDPack.h>
#include <blaze/math/shims/Clamp.h>
#include <blaze/math/simd/Max.h>
#include <blaze/math/simd/Min.h>
#include <blaze/math/simd/Set.h>
#include <blaze/math/typetraits/HasSIMDMax.h>
#include <blaze/math/typetraits/HasSIMDMin.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for the clamp() function.
// \ingroup functors
*/
template< typename DT >  // Type of the delimiters
struct Clamp
{
 public:
   //**********************************************************************************************
   /*!\brief Constructor of the Clamp functor.
   //
   // \param min The lower limit of the range.
   // \param max The upper limit of the range.
   */
   explicit inline Clamp( const DT& min, const DT& max )
      : min_( min )  // The lower delimiter
      , max_( max )  // The upper delimiter
   {}
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the clamp() function for the given object/value.
   //
   // \param a The given object/value.
   // \return The result of the clamp() function for the given object/value.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( const T& a ) const
   {
      return clamp( a, min_, max_ );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data type \a T.
   //
   // \return \a true in case SIMD is enabled for the data type \a T, \a false if not.
   */
   template< typename T >
   static constexpr bool simdEnabled() { return HasSIMDMin_v<T,DT> && HasSIMDMax_v<T,DT>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operation supports padding, i.e. whether it can deal with zeros.
   //
   // \return \a true in case padding is supported, \a false if not.
   */
   static constexpr bool paddingEnabled() { return true; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the clamp() function for the given SIMD vector.
   //
   // \param a The given SIMD vector.
   // \return The result of the clamp() function for the given SIMD vector.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const T& a ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK( T );
      return max( min( a, set( max_ ) ), set( min_ ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   DT min_;  //!< The lower delimiter.
   DT max_;  //!< The upper delimiter.
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
