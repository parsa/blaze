//=================================================================================================
/*!
//  \file blaze/util/algorithms/Max.h
//  \brief Headerfile for the generic max algorithm
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

#ifndef _BLAZE_UTIL_ALGORITHMS_MAX_H_
#define _BLAZE_UTIL_ALGORITHMS_MAX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/DisableMax.h>
#include <blaze/util/typetraits/IsSigned.h>
#include <blaze/util/typetraits/IsUnsigned.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  MAX ALGORITHMS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Maximum function for two values/objects.
// \ingroup algorithms
//
// \param a The first value/object.
// \param b The second value/object.
// \return The maximum of the two values/objects.
//
// This function determines the maximum of the two given values/objects by means of a less-than
// comparison. In case both values are equal, \a a is returned.
//
// The return type of the function is determined by the data types of the given arguments. In case
// two lvalues are given the function returns an lvalue reference. In case two rvalues are given
// it returns an rvalue reference. Else it returns by value.
*/
template< typename T1, typename T2
        , typename = EnableIf_t< !DisableMax_v< RemoveReference_t<T1> > &&
                                 !DisableMax_v< RemoveReference_t<T2> > &&
                                 !( IsSigned_v< RemoveReference_t<T1> > && IsUnsigned_v< RemoveReference_t<T2> > ) &&
                                 !( IsUnsigned_v< RemoveReference_t<T1> > && IsSigned_v< RemoveReference_t<T2> > ) > >
inline constexpr decltype(auto)
   max( T1&& a, T2&& b ) noexcept
{
   using std::forward;

   return ( a < b ) ? forward<T2>( b ) : forward<T1>( a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Maximum function for three values/objects.
// \ingroup algorithms
//
// \param a The first value/object.
// \param b The second value/object.
// \param c The third value/object.
// \return The maximum of the given values/objects.
//
// This function returns the maximum of the given data values/objects.
//
// The return type of the function is determined by the data types of the given arguments. In
// case three lvalues are given the function returns an lvalue reference. In case three rvalues
// are given it returns an rvalue reference. Else it returns by value.
*/
template< typename T1, typename T2, typename T3 >
inline constexpr decltype(auto)
   max( T1&& a, T2&& b, T3&& c ) noexcept
{
   using std::forward;
   using blaze::max;

   return max( max( forward<T1>( a ), forward<T2>( b ) ), forward<T3>( c ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Maximum function for at least four values/objects.
// \ingroup algorithms
//
// \param a The first value/object.
// \param b The second value/object.
// \param c The third value/object.
// \param args The pack of additional values/objects.
// \return The maximum of the given values/objects.
//
// This function returns the maximum of the given data values/objects.
//
// The return type of the function is determined by the data types of the given arguments. In
// case only lvalues are given the function returns an lvalue reference. In case only rvalues
// are given it returns an rvalue reference. Else it returns by value.
*/
template< typename T1, typename T2, typename T3, typename... Ts >
inline constexpr decltype(auto)
   max( T1&& a, T2&& b, T3&& c, Ts&&... args ) noexcept
{
   using std::forward;
   using blaze::max;

   return max( max( max( forward<T1>( a ), forward<T2>( b ) ), forward<T3>( c ) ), forward<Ts>( args )... );
}
//*************************************************************************************************

} // namespace blaze

#endif
