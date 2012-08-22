//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsObject.h
//  \brief Header file for the IsObject type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISOBJECT_H_
#define _BLAZE_UTIL_TYPETRAITS_ISOBJECT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_object.hpp>
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
/*!\brief Auxiliary helper struct for the IsObject type trait.
// \ingroup type_traits
*/
template< typename T >
struct IsObjectHelper
{
   //**********************************************************************************************
   enum { value = boost::is_object<T>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// This class tests whether the given template parameter \a T is an object type. All types are
// considered object types except references, \a void, and function types. If \a T is an object
// type, the \a value member enumeration is set to 1, the nested type definition \a Type is
// \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set to 0, \a Type
// is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsObject<int>::value                   // Evaluates to 1
   blaze::IsObject<int*>::Type                   // Results in TrueType
   blaze::IsObject<int (*)(void)>                // Is derived from TrueType
   blaze::IsObject<int (MyClass::*)(void)const>  // Also derived from TrueType
   blaze::IsObject<int&>::value                  // Evaluates to 0
   blaze::IsObject<const void>::Type             // Results in FalseType
   blaze::IsObject<int (double)>                 // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsObject : public IsObjectHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsObjectHelper<T>::value };
   typedef typename IsObjectHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
