//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsPointer.h
//  \brief Header file for the IsPointer type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISPOINTER_H_
#define _BLAZE_UTIL_TYPETRAITS_ISPOINTER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_pointer.hpp>
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
/*!\brief Auxiliary helper struct for the IsPointer type trait.
// \ingroup type_traits
*/
template< typename T >
struct IsPointerHelper
{
   //**********************************************************************************************
   enum { value = boost::is_pointer<T>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// This class tests whether the given template parameter is a pointer type (including function
// pointers, but excluding pointers to members) or not. If it is a pointer type, the \a value
// member enumeration is set to 1, the nested type definition \a Type is \a TrueType, and the
// class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.

   \code
   blaze::IsPointer<char* const>::value      // Evaluates to 1
   blaze::IsPointer<volatile float*>::Type   // Results in TrueType
   blaze::IsPointer<int (*)(long)>           // Is derived from TrueType
   blaze::IsPointer<int>::value              // Evaluates to 0
   blaze::IsPointer<int MyClass::*>::Type    // Results in FalseType
   blaze::IsPointer<int (MyClass::*)(long)>  // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsPointer : public IsPointerHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsPointerHelper<T>::value };
   typedef typename IsPointerHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
