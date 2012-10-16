//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsTemporary.h
//  \brief Header file for the IsTemporary type trait class
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISTEMPORARY_H_
#define _BLAZE_MATH_TYPETRAITS_ISTEMPORARY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/IsReference.h>
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
/*!\brief Auxiliary helper struct for the IsTemporary type trait.
// \ingroup math_type_traits
*/
template< typename T >
struct IsTemporaryHelper
{
   //**********************************************************************************************
   enum { value = !IsReference<T>::value && !IsNumeric<T>::value && !IsExpression<T>::value };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check whether the given type is a temporary vector or matrix type.
// \ingroup math_type_traits
//
// This type trait class tests whether the given type is a temporary vector or matrix type,
// i.e. can be used for a temporary vector or matrix. In case the given type can be used as
// temporary, the \a value member enumeration is set to 1, the nested type definition \a Type
// is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set to 0,
// \a Type is \a FalseType, and the class derives from \a FalseType.
*/
template< typename T >
struct IsTemporary : public IsTemporaryHelper<T>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsTemporaryHelper<T>::value };
   typedef typename IsTemporaryHelper<T>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
