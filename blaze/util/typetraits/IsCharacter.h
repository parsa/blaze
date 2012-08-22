//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsCharacter.h
//  \brief Header file for the IsCharacter type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISCHARACTER_H_
#define _BLAZE_UTIL_TYPETRAITS_ISCHARACTER_H_


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
/*!\brief Compile time check for character types.
// \ingroup type_traits
//
// This type trait tests whether or not the given template parameter is a character type
// (i.e., either char, signed char, unsigned char, or wchar_t, possibly cv-qualified). In
// case the type is a character type (ignoring the cv-qualifiers), the \a value member
// enumeration is set to 1, the nested type definition \a Type is \a TrueType, and the
// class derives from \a TrueType. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.

   \code
   blaze::IsCharacter<char>::value                // Evaluates to 1
   blaze::IsCharacter<const unsigned char>::Type  // Results in TrueType
   blaze::IsCharacter<const volatile wchar_t>     // Is derived from TrueType
   blaze::IsCharacter<unsigned short>::value      // Evaluates to 0
   blaze::IsCharacter<const int>::Type            // Results in FalseType
   blaze::IsCharacter<volatile long>              // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsCharacter : public FalseType
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
//! Specialization of the IsCharacter type trait for the plain 'char' type.
template<>
struct IsCharacter<char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'const char'.
template<>
struct IsCharacter<const char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'volatile char'.
template<>
struct IsCharacter<volatile char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'const volatile char'.
template<>
struct IsCharacter<const volatile char> : public TrueType
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
//! Specialization of the IsCharacter type trait for the plain 'signed char' type.
template<>
struct IsCharacter<signed char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'const signed char'.
template<>
struct IsCharacter<const signed char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'volatile signed char'.
template<>
struct IsCharacter<volatile signed char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'const volatile signed char'.
template<>
struct IsCharacter<const volatile signed char> : public TrueType
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
//! Specialization of the IsCharacter type trait for the plain 'unsigned char' type.
template<>
struct IsCharacter<unsigned char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'const unsigned char'.
template<>
struct IsCharacter<const unsigned char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'volatile unsigned char'.
template<>
struct IsCharacter<volatile unsigned char> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'const volatile unsigned char'.
template<>
struct IsCharacter<const volatile unsigned char> : public TrueType
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
//! Specialization of the IsCharacter type trait for the plain 'wchar_t' type.
template<>
struct IsCharacter<wchar_t> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'const wchar_t'.
template<>
struct IsCharacter<const wchar_t> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'volatile wchar_t'.
template<>
struct IsCharacter<volatile wchar_t> : public TrueType
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
//! Specialization of the IsCharacter type trait for 'const volatile wchar_t'.
template<>
struct IsCharacter<const volatile wchar_t> : public TrueType
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
