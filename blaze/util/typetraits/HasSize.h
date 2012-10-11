//=================================================================================================
/*!
//  \file blaze/util/typetraits/HasSize.h
//  \brief Header file for the HasSize type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_HASSIZE_H_
#define _BLAZE_UTIL_TYPETRAITS_HASSIZE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/FalseType.h>
#include <blaze/util/SelectType.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS HASSIZE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time size check.
// \ingroup type_traits
//
// This class offers the possibility to test the size of a type at compile time. If the type
// \a T is exactly \a Size bytes large, the \a value member enumeration is set to 1, the nested
// type definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to 0, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::HasSize<int,4>::value              // Evaluates to 1 (on most architectures)
   blaze::HasSize<float,4>::Type             // Results in TrueType (on most architectures)
   blaze::HasSize<const double,8>            // Is derived from TrueType (on most architectures)
   blaze::HasSize<volatile double,2>::value  // Evaluates to 0
   blaze::HasSize<const char,8>::Type        // Results in FalseType
   blaze::HasSize<unsigned char,4>           // Is derived from FalseType
   \endcode
*/
template< typename T, size_t Size >
struct HasSize : public SelectType< sizeof( T ) == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = ( sizeof( T ) == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class ia a partial specialization of the HasSize template for the type \a void. This
// specialization assumes that an object of type \a void has a size of 0. Therefore \a value
// is set to 1, \a Type is \a TrueType, and the class derives from \a TrueType only if the
// \a Size template argument is 0. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< size_t Size >
struct HasSize<void,Size> : public SelectType< 0 == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( 0 == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class ia a partial specialization of the HasSize template for constant \a void. This
// specialization assumes that an object of type \a void has a size of 0. Therefore \a value
// is set to 1, \a Type is \a TrueType, and the class derives from \a TrueType only if the
// \a Size template argument is 0. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< size_t Size >
struct HasSize<const void,Size> : public SelectType< 0 == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( 0 == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class ia a partial specialization of the HasSize template for volatile \a void. This
// specialization assumes that an object of type \a void has a size of 0. Therefore \a value
// is set to 1, \a Type is \a TrueType, and the class derives from \a TrueType only if the
// \a Size template argument is 0. Otherwise \a value is set to 0, \a Type is \a FalseType,
// and the class derives from \a FalseType.
*/
template< size_t Size >
struct HasSize<volatile void,Size> : public SelectType< 0 == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( 0 == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Partial specialization of the compile time size constraint.
// \ingroup type_traits
//
// This class ia a partial specialization of the HasSize template for constant volatile \a void.
// This specialization assumes that an object of type \a void has a size of 0. Therefore \a value
// is set to 1, \a Type is \a TrueType, and the class derives from \a TrueType only if the \a Size
// template argument is 0. Otherwise \a value is set to 0, \a Type is \a FalseType, and the class
// derives from \a FalseType.
*/
template< size_t Size >
struct HasSize<const volatile void,Size> : public SelectType< 0 == Size, TrueType, FalseType >::Type
{
 public:
   //**********************************************************************************************
   enum { value = ( 0 == Size ) };
   typedef typename SelectType<value,TrueType,FalseType>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS HAS1BYTE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time size check.
// \ingroup type_traits
//
// This type trait offers the possibility to test whether a given type has a size of exactly
// one byte. If the type \a T has one byte, the \a value member enumeration is set to 1, the
// nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives from
// \a FalseType.

   \code
   blaze::Has1Byte<const char>::value       // Evaluates to 1 (on most architectures)
   blaze::Has1Byte<unsigned char>::Type     // Results in TrueType (on most architectures)
   blaze::Has1Byte<signed char>             // Is derived from TrueType (on most architectures)
   blaze::Has1Byte<volatile double>::value  // Evaluates to 0
   blaze::Has1Byte<const float>::Type       // Results in FalseType
   blaze::Has1Byte<unsigned short>          // Is derived from FalseType
   \endcode
*/
template< typename T >
struct Has1Byte : public HasSize<T,1UL>
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS HAS2BYTES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time size check.
// \ingroup type_traits
//
// This type trait offers the possibility to test whether a given type has a size of exactly
// two bytes. If the type \a T has two bytes, the \a value member enumeration is set to 1, the
// nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives from
// \a FalseType.

   \code
   blaze::Has2Bytes<const short>::value      // Evaluates to 1 (on most architectures)
   blaze::Has2Bytes<unsigned short>::Type    // Results in TrueType (on most architectures)
   blaze::Has2Bytes<volatile short>          // Is derived from TrueType (on most architectures)
   blaze::Has2Bytes<volatile double>::value  // Evaluates to 0
   blaze::Has2Bytes<const float>::Type       // Results in FalseType
   blaze::Has2Bytes<unsigned int>            // Is derived from FalseType
   \endcode
*/
template< typename T >
struct Has2Bytes : public HasSize<T,2UL>
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS HAS4BYTES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time size check.
// \ingroup type_traits
//
// This type trait offers the possibility to test whether a given type has a size of exactly
// four bytes. If the type \a T has four bytes, the \a value member enumeration is set to 1,
// the nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives from
// \a FalseType.

   \code
   blaze::Has4Bytes<const int>::value        // Evaluates to 1 (on most architectures)
   blaze::Has4Bytes<unsigned int>::Type      // Results in TrueType (on most architectures)
   blaze::Has4Bytes<volatile float>          // Is derived from TrueType (on most architectures)
   blaze::Has4Bytes<volatile double>::value  // Evaluates to 0
   blaze::Has4Bytes<const float>::Type       // Results in FalseType
   blaze::Has4Bytes<short>                   // Is derived from FalseType
   \endcode
*/
template< typename T >
struct Has4Bytes : public HasSize<T,4UL>
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS HAS8BYTES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time size check.
// \ingroup type_traits
//
// This type trait offers the possibility to test whether a given type has a size of exactly
// four bytes. If the type \a T has four bytes, the \a value member enumeration is set to 1,
// the nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives from
// \a FalseType.

   \code
   blaze::Has8Bytes<double>::value        // Evaluates to 1 (on most architectures)
   blaze::Has8Bytes<const double>::Type   // Results in TrueType (on most architectures)
   blaze::Has8Bytes<volatile double>      // Is derived from TrueType (on most architectures)
   blaze::Has8Bytes<unsigned int>::value  // Evaluates to 0
   blaze::Has8Bytes<const float>::Type    // Results in FalseType
   blaze::Has8Bytes<volatile short>       // Is derived from FalseType
   \endcode
*/
template< typename T >
struct Has8Bytes : public HasSize<T,8UL>
{};
//*************************************************************************************************

} // namespace blaze

#endif
