//=================================================================================================
/*!
//  \file blaze/util/typetraits/Extent.h
//  \brief Header file for the Extent type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_EXTENT_H_
#define _BLAZE_UTIL_TYPETRAITS_EXTENT_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for the size of array bounds.
// \ingroup type_traits
//
// Via this type trait it is possible to query at compile time for the size of a particular
// array extent. In case the given template argument is an array type with a rank greater
// than N, the \a value member enumeration is set to the number of elements of the N'th
// array dimension. In all other cases, and especially in case the N'th array dimension
// is incomplete, \a value is set to 0.

   \code
   blaze::Extent< int[4], 0 >::value            // Evaluates to 4
   blaze::Extent< int[2][3][4], 0 >::value      // Evaluates to 2
   blaze::Extent< int[2][3][4], 1 >::value      // Evaluates to 3
   blaze::Extent< int[2][3][4], 2 >::value      // Evaluates to 4
   blaze::Extent< int[][2], 0 >::value          // Evaluates to 0
   blaze::Extent< int[][2], 1 >::value          // Evaluates to 2
   blaze::Extent< int*, 0 >::value              // Evaluates to 0
   blaze::Extent< std::vector<int>, 0 >::value  // Evaluates to 0 (std::vector is NOT an array type)
   \endcode
*/
template< typename T, unsigned int N >
struct Extent
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = 0 };
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Partial specialization of the Extent type trait for empty array extents.
template< typename T, unsigned int N >
struct Extent<T[],N>
{
 public:
   //**********************************************************************************************
   enum { value = Extent<T>::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Partial specialization of the Extent type trait for non-empty array extents.
template< typename T, unsigned int N, unsigned int E >
struct Extent<T[E],N>
{
 public:
   //**********************************************************************************************
   enum { value = Extent<T>::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Terminating partial specialization of the Extent type trait for empty array extents.
template< typename T >
struct Extent<T[],0UL>
{
 public:
   //**********************************************************************************************
   enum { value = 0 };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Terminating partial specialization of the Extent type trait for non-empty array extents.
template< typename T, unsigned int E >
struct Extent<T[E],0U>
{
 public:
   //**********************************************************************************************
   enum { value = E };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
