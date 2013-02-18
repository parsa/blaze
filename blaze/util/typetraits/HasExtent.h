//=================================================================================================
/*!
//  \file blaze/util/typetraits/HasExtent.h
//  \brief Header file for the HasExtent type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_HASEXTENT_H_
#define _BLAZE_UTIL_TYPETRAITS_HASEXTENT_H_


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
/*!\brief Compile time check for array extents.
// \ingroup type_traits
//
// This type trait tests whether or not the given template argument has any array extents
// (empty or non-empty) and determines the number of array extents. In case the type has any
// array extent, the \a value member enumeration is set to the total number of array extents,
// the nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to 0, \a Type is \a FalseType, and the class derives from
// \a FalseType.

   \code
   blaze::HasExtent<int[3]>::value      // Evaluates to 1
   blaze::HasExtent<const int[]>::Type  // Results in TrueType
   blaze::HasExtent<int[][3]>           // Is derived from TrueType
   blaze::HasExtent<int>::value         // Evaluates to 0
   blaze::HasExtent<int const*>::Type   // Results in FalseType
   blaze::HasExtent<int volatile**>     // Is derived from FalseType
   \endcode
*/
template< typename T >
struct HasExtent : public FalseType
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
//! Specialization of the HasExtent type trait for empty array extents.
template< typename T >
struct HasExtent<T[]> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 + HasExtent<T>::value };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the HasExtent type trait for non-empty array extents.
template< typename T, unsigned int N >
struct HasExtent<T[N]> : public TrueType
{
 public:
   //**********************************************************************************************
   enum { value = 1 + HasExtent<T>::value };
   typedef TrueType  Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
