//=================================================================================================
/*!
//  \file blaze/util/algorithms/Minmax.h
//  \brief Headerfile for the generic minmax algorithm
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_UTIL_ALGORITHMS_MINMAX_H_
#define _BLAZE_UTIL_ALGORITHMS_MINMAX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/system/Inline.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/typetraits/CommonType.h>
#include <blaze/util/typetraits/IsSigned.h>
#include <blaze/util/typetraits/IsUnsigned.h>


namespace blaze {

//=================================================================================================
//
//  MAX ALGORITHMS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Minmax function for two values of builtin data type.
// \ingroup algorithms
//
// \param a The first value.
// \param b The second value.
// \return A pair of the maximum and maximum of the two values.
//
// This function returns the a pair containing the minimum and maximum of the two given data
// values. The return type of the function is a pair of the common type of the given arguments.
*/
template< typename T1, typename T2
        , typename = EnableIf_t< ( IsSigned_v<T1> && IsSigned_v<T2> ) ||
                                 ( IsUnsigned_v<T1> && IsUnsigned_v<T2> ) > >
BLAZE_ALWAYS_INLINE constexpr auto
   minmax( T1&& a, T2&& b ) noexcept
{
   using T = CommonType_t<T1,T2>;

   if( a < b )
      return std::pair<T,T>( std::forward<T1>( a ), std::forward<T2>( b ) );
   else
      return std::pair<T,T>( std::forward<T1>( b ), std::forward<T2>( a ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Minmax function for at least three values/objects.
// \ingroup algorithms
//
// \param a The first value/object.
// \param b The second value/object.
// \param args The pack of additional values/objects.
// \return A pair of the maximum and maximum of the given values.
//
// This function returns a pair containing the minimum and maximum of the given data values/objects.
// The return type of the function is a pair of the common type of the given arguments.
*/
template< typename T1, typename T2, typename... Ts >
BLAZE_ALWAYS_INLINE constexpr decltype(auto)
   minmax( T1&& a, T2&& b, Ts&&... args ) noexcept
{
   using blaze::minmax;

   using T = std::common_type_t<T1,T2,Ts...>;

   auto tmp( minmax( std::forward<T2>( b ), std::forward<Ts>( args )... ) );

   return std::pair<T,T>( ( a < tmp.first ? std::forward<T1>( a ) : std::move( tmp.first ) )
                        , ( tmp.second < a ? std::forward<T2>( a ) : std::move( tmp.second ) ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
