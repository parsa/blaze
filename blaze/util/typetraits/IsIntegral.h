//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsIntegral.h
//  \brief Header file for the IsIntegral type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISINTEGRAL_H_
#define _BLAZE_UTIL_TYPETRAITS_ISINTEGRAL_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_integral.hpp>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsIntegral type trait.
// \ingroup type_traits
*/
template< typename T >
struct IsIntegralHelper
{
   //**********************************************************************************************
   enum { value = boost::is_integral<T>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for integral data types.
// \ingroup type_traits
//
// This type trait tests whether or not the given template parameter is an integral data
// type. In case the type is an integral data type, the \a value member enumeration is set
// to 1, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class
// derives from \a FalseType.

   \code
   blaze::IsIntegral<int>::value            // Evaluates to 1
   blaze::IsIntegral<const char>::Type      // Results in TrueType (char is an integral data type)
   blaze::IsIntegral<volatile short>        // Is derived from TrueType
   blaze::IsIntegral<float>::value          // Evaluates to 0
   blaze::IsIntegral<const double>::Type    // Results in FalseType
   blaze::IsIntegral<volatile long double>  // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsIntegral : public IsIntegralHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsIntegralHelper<T>::value };
   typedef typename IsIntegralHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
