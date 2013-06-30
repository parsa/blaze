//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsBaseOf.h
//  \brief Header file for the IsBaseOf type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISBASEOF_H_
#define _BLAZE_UTIL_TYPETRAITS_ISBASEOF_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/is_base_of.hpp>
#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/TypeTraits/RemoveCV.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsBaseOf type trait.
// \ingroup type_traits
*/
template< typename Base, typename Derived >
struct IsBaseOfHelper
{
   //**********************************************************************************************
   enum { value = boost::is_base_of<typename RemoveCV<Base>::Type,
                                    typename RemoveCV<Derived>::Type>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time analysis of an inheritance relationship.
// \ingroup type_traits
//
// This type trait tests for an inheritance relationship between the two types \a Base and
// \a Derived. If \a Derived is a type derived from \a Base or the same type as \a Base the
// \a value member enumeration is set to 1, the nested type definition \a Type is \a TrueType,
// and the class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is
// \a FalseType, and the class derives from \a FalseType.

   \code
   class A { ... };
   class B : public A { ... };
   class C { ... };

   blaze::IsBaseOf<A,B>::value  // Evaluates to 1
   blaze::IsBaseOf<A,B>::Type   // Results in TrueType
   blaze::IsBaseOf<A,B>         // Is derived from TrueType
   blaze::IsBaseOf<A,C>::value  // Evaluates to 0
   blaze::IsBaseOf<B,A>::Type   // Results in FalseType
   blaze::IsBaseOf<B,A>         // Is derived from FalseType
   \endcode
*/
template< typename Base, typename Derived >
class IsBaseOf : public IsBaseOfHelper<Base,Derived>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsBaseOfHelper<Base,Derived>::value };
   typedef typename IsBaseOfHelper<Base,Derived>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
