//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsEmpty.h
//  \brief Header file for the IsEmpty type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISEMPTY_H_
#define _BLAZE_UTIL_TYPETRAITS_ISEMPTY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_empty.hpp>
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
/*!\brief Auxiliary helper struct for the IsEmpty type trait.
// \ingroup type_traits
*/
template< typename T >
struct IsEmptyHelper
{
   //**********************************************************************************************
   enum { value = boost::is_empty<T>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time type check.
// \ingroup type_traits
//
// This class tests whether the given template parameter is an empty class type, i.e. a type
// without member data and virtual functions. If it is an empty class type, the \a value
// member enumeration is set to 1, the nested type definition \a Type is \a TrueType, and the
// class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.

   \code
   class A {};
   class B { int i; };

   blaze::IsEmpty<A>::value           // Evaluates to 1
   blaze::IsEmpty<A volatile>::Type   // Results in TrueType
   blaze::IsEmpty<A const>            // Is derived from TrueType
   blaze::IsEmpty<int>::value         // Evaluates to 0
   blaze::IsEmpty<std::string>::Type  // Results in FalseType
   blaze::IsEmpty<B>                  // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsEmpty : public IsEmptyHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsEmptyHelper<T>::value };
   typedef typename IsEmptyHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
