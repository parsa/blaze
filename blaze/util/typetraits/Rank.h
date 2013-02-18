//=================================================================================================
/*!
//  \file blaze/util/typetraits/Rank.h
//  \brief Header file for the Rank type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_RANK_H_
#define _BLAZE_UTIL_TYPETRAITS_RANK_H_


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
/*!\brief Compile time check for array ranks.
// \ingroup type_traits
//
// This type trait determines the rank of the given template argument. In case the given type
// is an array type, the nested \a value member enumeration is set to the number of dimensions
// of \a T. Otherwise \a value is set to 0.

   \code
   blaze::Rank< int[] >::value               // Evaluates to 1
   blaze::Rank< int[3] >::value              // Evaluates to 1
   blaze::Rank< const int[2][3][4] >::value  // Evaluates to 3
   blaze::Rank< int[][3] >::value            // Evaluates to 2
   blaze::Rank< int const* >::value          // Evaluates to 0
   blaze::Rank< std::vector<int> >::value    // Evaluates to 0
   \endcode
*/
template< typename T >
struct Rank
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
//! Specialization of the Rank type trait for empty arrays.
template< typename T >
struct Rank<T[]>
{
 public:
   //**********************************************************************************************
   enum { value = 1 + Rank<T>::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
//! Specialization of the Rank type trait for non-empty arrays.
template< typename T, unsigned int N >
struct Rank<T[N]>
{
 public:
   //**********************************************************************************************
   enum { value = 1 + Rank<T>::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
