//=================================================================================================
/*!
//  \file blaze/util/typetraits/RemoveConst.h
//  \brief Header file for the RemoveConst type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_REMOVECONST_H_
#define _BLAZE_UTIL_TYPETRAITS_REMOVECONST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_const.hpp>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Removal of const-qualifiers.
// \ingroup type_traits
//
// The RemoveConst type trait removes all top level 'const' qualifiers from the given type \a T.

   \code
   blaze::RemoveConst<short>::Type               // Results in 'short'
   blaze::RemoveConst<const double>::Type        // Results in 'double'
   blaze::RemoveConst<const volatile int>::Type  // Results in 'volatile int'
   blaze::RemoveConst<int const*>::Type          // Results in 'const int*'
   blaze::RemoveConst<int const* const>::Type    // Results in 'const int*'
   blaze::RemoveConst<int const&>::Type          // Results in 'const int&'
   \endcode
*/
template< typename T >
struct RemoveConst
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename boost::remove_const<T>::type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
