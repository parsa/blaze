//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsSame.h
//  \brief Header file for the IsSame and IsStrictlySame type traits
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISSAME_H_
#define _BLAZE_UTIL_TYPETRAITS_ISSAME_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_cv.hpp>
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
/*!\brief Compile time type relationship analysis.
// \ingroup type_traits
//
// This class tests if the two data types \a A and \a B are equal. For this type comparison,
// the cv-qualifiers of both data types are not ignored. If \a A and \a B are the same data
// type, then the \a value member enumeration is set to 1, the nested type definition \a Type
// is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set to 0,
// \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsStrictlySame<int,int>::value                   // Evaluates to 1
   blaze::IsStrictlySame<const double,const double>::Type  // Results in TrueType
   blaze::IsStrictlySame<volatile float,volatile float>    // Is derived from TrueType
   blaze::IsStrictlySame<char,wchar_t>::value              // Evaluates to 0
   blaze::IsStrictlySame<int,const int>::Type              // Results in FalseType
   blaze::IsStrictlySame<float,volatile float>             // Is derived from FalseType
   \endcode
*/
template< typename A, typename B >
struct IsStrictlySame : public FalseType
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
//! Specialization of the IsStrictlySame class template for a single, matching data type.
template< typename T >
struct IsStrictlySame<T,T> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsSame type trait.
// \ingroup type_traits
*/
template< typename A, typename B >
struct IsSameHelper
{
 private:
   //**********************************************************************************************
   typedef typename boost::remove_cv<A>::type  T1;
   typedef typename boost::remove_cv<B>::type  T2;
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   enum { value = IsStrictlySame<T1,T2>::value };
   typedef typename IsStrictlySame<T1,T2>::Type  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Type relationship analysis.
// \ingroup type_traits
//
// This class tests if the two data types \a A and \a B are equal. For this type comparison,
// the cv-qualifiers of both data types are ignored. If \a A and \a B are the same data type
// (ignoring the cv-qualifiers), then the \a value member enumeration is set to 1, the nested
// type definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to 0, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsSame<int,int>::value               // Evaluates to 1
   blaze::IsSame<int,const int>::Type          // Results in TrueType
   blaze::IsSame<float,volatile float>         // Is derived from TrueType
   blaze::IsSame<char,wchar_t>::value          // Evaluates to 0
   blaze::IsSame<char,volatile float>::Type    // Results in FalseType
   blaze::IsSame<int,double>                   // Is derived from FalseType
   \endcode
*/
template< typename A, typename B >
struct IsSame : public IsSameHelper<A,B>::Type
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = IsSameHelper<A,B>::value };
   typedef typename IsSameHelper<A,B>::Type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
