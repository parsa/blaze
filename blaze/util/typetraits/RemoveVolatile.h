//=================================================================================================
/*!
//  \file blaze/util/typetraits/RemoveVolatile.h
//  \brief Header file for the RemoveVolatile type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_REMOVEVOLATILE_H_
#define _BLAZE_UTIL_TYPETRAITS_REMOVEVOLATILE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_volatile.hpp>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Removal of volatile-qualifiers.
// \ingroup type_traits
//
// The RemoveVolatile type trait removes all top level 'volatile' qualifiers from the given
// type \a T.

   \code
   blaze::RemoveVolatile<short>::Type                   // Results in 'short'
   blaze::RemoveVolatile<volatile double>::Type         // Results in 'double'
   blaze::RemoveVolatile<const volatile int>::Type      // Results in 'const int'
   blaze::RemoveVolatile<int volatile*>::Type           // Results in 'const int*'
   blaze::RemoveVolatile<int volatile* volatile>::Type  // Results in 'const int*'
   blaze::RemoveVolatile<int volatile&>::Type           // Results in 'const int&'
   \endcode
*/
template< typename T >
struct RemoveVolatile
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename boost::remove_volatile<T>::type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
