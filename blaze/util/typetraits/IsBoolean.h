//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsBoolean.h
//  \brief Header file for the IsBoolean type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISBOOLEAN_H_
#define _BLAZE_UTIL_TYPETRAITS_ISBOOLEAN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/FalseType.h>
#include <blaze/util/TrueType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for boolean types.
// \ingroup type_traits
//
// This type trait tests whether or not the given template parameter is of boolean type. In
// case the type is a boolean (ignoring the cv-qualifiers), the \a value member enumeration
// is set to 1, the nested type definition \a Type is \a TrueType, and the class derives
// from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class
// derives from \a FalseType.

   \code
   blaze::IsBoolean<bool>::value          // Evaluates to 1
   blaze::IsBoolean<const bool>::Type     // Results in TrueType
   blaze::IsBoolean<const volatile bool>  // Is derived from TrueType
   blaze::IsBoolean<float>::value         // Evaluates to 0 (float is not a boolean)
   blaze::IsBoolean<const int>::Type      // Results in FalseType
   blaze::IsBoolean<volatile short>       // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsBoolean : public FalseType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = 0 };
   typedef FalseType  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsBoolean type trait for the plain 'bool' type.
template<>
struct IsBoolean<bool> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsBoolean type trait for 'const bool'.
template<>
struct IsBoolean<const bool> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsBoolean type trait for 'volatile bool'.
template<>
struct IsBoolean<volatile bool> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the IsBoolean type trait for 'const volatile bool'
template<>
struct IsBoolean<const volatile bool> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
