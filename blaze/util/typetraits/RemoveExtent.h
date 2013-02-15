//=================================================================================================
/*!
//  \file blaze/util/typetraits/RemoveExtent.h
//  \brief Header file for the RemoveExtent type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_REMOVEEXTENT_H_
#define _BLAZE_UTIL_TYPETRAITS_REMOVEEXTENT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_extent.hpp>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Removal of the top level array extent.
// \ingroup type_traits
//
// The RemoveExtent type trait removes the top level array extent from the given type \a T.

   \code
   blaze::RemoveExtent<int>::Type           // Results in 'int'
   blaze::RemoveExtent<int const[2]>::Type  // Results in 'int const'
   blaze::RemoveExtent<int[2][4]>::Type     // Results in 'int[4]'
   blaze::RemoveExtent<int[][2]>::Type      // Results in 'int[2]'
   blaze::RemoveExtent<int const*>::Type    // Results in 'int const*'
   \endcode
*/
template< typename T >
struct RemoveExtent
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename boost::remove_extent<T>::type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
