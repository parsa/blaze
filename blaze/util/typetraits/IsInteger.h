//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsInteger.h
//  \brief Header file for the IsInteger type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISINTEGER_H_
#define _BLAZE_UTIL_TYPETRAITS_ISINTEGER_H_


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
/*!\brief Compile time check for integer types.
// \ingroup type_traits
//
// This type trait tests whether or not the given template parameter is an integer type (i.e.,
// either (signed) int or unsigned int, possibly cv-qualified). In case the type is an integer
// type (ignoring the cv-qualifiers), the \a value member enumeration is set to 1, the nested
// type definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to 0, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsInteger<int>::value                 // Evaluates to 1
   blaze::IsInteger<const unsigned int>::Type   // Results in TrueType
   blaze::IsInteger<const volatile signed int>  // Is derived from TrueType
   blaze::IsInteger<unsigned short>::value      // Evaluates to 0
   blaze::IsInteger<const long>::Type           // Results in FalseType
   blaze::IsInteger<volatile float>             // Is derived from FalseType
   \endcode

// Note the difference between the IsInteger and IsIntegral type traits: Whereas the IsInteger
// type trait specifically tests whether the given data type is either int or unsigned int
// (possibly cv-qualified), the IsIntegral type trait tests whether the given template argument
// is an integral data type (char, short, int, long, etc.).
*/
template< typename T >
struct IsInteger : public FalseType
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
//! Specialization of the IsInteger type trait for the plain 'int' type.
template<>
struct IsInteger<int> : public TrueType
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
//! Specialization of the IsInteger type trait for 'const int'.
template<>
struct IsInteger<const int> : public TrueType
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
//! Specialization of the IsInteger type trait for 'volatile int'.
template<>
struct IsInteger<volatile int> : public TrueType
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
//! Specialization of the IsInteger type trait for 'const volatile int'.
template<>
struct IsInteger<const volatile int> : public TrueType
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
//! Specialization of the IsInteger type trait for the plain 'unsigned int' type.
template<>
struct IsInteger<unsigned int> : public TrueType
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
//! Specialization of the IsInteger type trait for 'const unsigned int'.
template<>
struct IsInteger<const unsigned int> : public TrueType
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
//! Specialization of the IsInteger type trait for 'volatile unsigned int'.
template<>
struct IsInteger<volatile unsigned int> : public TrueType
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
//! Specialization of the IsInteger type trait for 'const volatile unsigned int'.
template<>
struct IsInteger<const volatile unsigned int> : public TrueType
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
