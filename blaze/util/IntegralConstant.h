//=================================================================================================
/*!
//  \file blaze/util/IntegralConstant.h
//  \brief Header file for the IntegralConstant class template
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

#ifndef _BLAZE_UTIL_INTEGRALCONSTANT_H_
#define _BLAZE_UTIL_INTEGRALCONSTANT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <type_traits>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for a compile time constant integral value.
// \ingroup util
//
// The IntegralConstant class template represents a generic wrapper for a compile time constant
// integral value. The value of an IntegralConstant can be accessed via the nested \a value (which
// is guaranteed to be of type \a T), the type can be accessed via the nested type definition
// \a ValueType.

   \code
   using namespace blaze;

   IntegralConstant<int,3>::value        // Evaluates to 3
   IntegralConstant<long,5L>::ValueType  // Results in long
   \endcode
*/
template< typename T, T N >
struct IntegralConstant
   : public std::integral_constant<T,N>
{
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using ValueType = T;
   using Type = IntegralConstant<T,N>;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Generic wrapper for a compile time constant boolean value.
// \ingroup util
//
// The BoolConstant alias template represents a generic wrapper for a compile time constant
// boolean value. The value of a BoolConstant can be accessed via the nested \a value (which
// is guaranteed to be of type \c bool), the type can be accessed via the nested type definition
// \a ValueType.

   \code
   using namespace blaze;

   BoolConstant<true>::value       // Evaluates to true
   BoolConstant<false>::ValueType  // Results in bool
   \endcode
*/
template< bool B >
using BoolConstant = IntegralConstant<bool,B>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time integral constant wrapper for \a bool.
// \ingroup mpl
//
// The Bool_t alias template represents an integral wrapper for a compile time constant
// expression of type \a bool. The value of a Bool_t can be accessed via the nested \a value
// (which is guaranteed to be of type \a bool), the type can be accessed via the nested type
// definition \a ValueType.

   \code
   using namespace blaze;

   Bool_t<true>::value       // Evaluates to true
   Bool_t<false>::ValueType  // Results in bool
   \endcode
*/
template< bool B >
using Bool_t = IntegralConstant<bool,B>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time integral constant wrapper for \a char.
// \ingroup mpl
//
// The Char_t alias template represents an integral wrapper for a compile time constant
// expression of type \a char. The value of an Char_t can be accessed via the nested \a value
// (which is guaranteed to be of type \a char), the type can be accessed via the nested type
// definition \a ValueType.

   \code
   using namespace blaze;

   Char_t<3>::value      // Evaluates to 3
   Char_t<5>::ValueType  // Results in char
   \endcode
*/
template< char N >
using Char_t = IntegralConstant<char,N>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time integral constant wrapper for \a int.
// \ingroup mpl
//
// The Int_t alias template represents an integral wrapper for a compile time constant
// expression of type \a int. The value of an Int_t can be accessed via the nested \a value
// (which is guaranteed to be of type \a int), the type can be accessed via the nested type
// definition \a ValueType.

   \code
   using namespace blaze;

   Int_t<3>::value      // Evaluates to 3
   Int_t<5>::ValueType  // Results in int
   \endcode
*/
template< int N >
using Int_t = IntegralConstant<int,N>;
//*************************************************************************************************

} // namespace blaze

#endif
