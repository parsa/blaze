//=================================================================================================
/*!
//  \file blaze/util/typetraits/RemoveAllExtents.h
//  \brief Header file for the RemoveAllExtents type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_REMOVEALLEXTENTS_H_
#define _BLAZE_UTIL_TYPETRAITS_REMOVEALLEXTENTS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_all_extents.hpp>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Removal of all array extents.
// \ingroup type_traits
//
// The RemoveAllExtents type trait removes all array extents from the given type \a T.

   \code
   blaze::RemoveAllExtents<int>::Type           // Results in 'int'
   blaze::RemoveAllExtents<int const[2]>::Type  // Results in 'int const'
   blaze::RemoveAllExtents<int[2][4]>::Type     // Results in 'int'
   blaze::RemoveAllExtents<int[][2]>::Type      // Results in 'int'
   blaze::RemoveAllExtents<int[2][3][4]>::Type  // Results in 'int'
   blaze::RemoveAllExtents<int const*>::Type    // Results in 'int const*'
   \endcode
*/
template< typename T >
struct RemoveAllExtents
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename boost::remove_all_extents<T>::type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
