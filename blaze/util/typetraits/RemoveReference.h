//=================================================================================================
/*!
//  \file blaze/util/typetraits/RemoveReference.h
//  \brief Header file for the RemoveReference type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_REMOVEREFERENCE_H_
#define _BLAZE_UTIL_TYPETRAITS_REMOVEREFERENCE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/remove_reference.hpp>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Removal of reference modifiers.
// \ingroup type_traits
//
// The RemoveCV type trait removes any reference modifiers from the given type \a T.

   \code
   blaze::RemoveCV<int>::Type             // Results in 'int'
   blaze::RemoveCV<const int&>::Type      // Results in 'const int'
   blaze::RemoveCV<volatile int&&>::Type  // Results in 'volatile int'
   blaze::RemoveCV<int*>::Type            // Results in 'int*'
   blaze::RemoveCV<int*&>::Type           // Results in 'int*'
   \endcode
*/
template< typename T >
struct RemoveReference
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename boost::remove_reference<T>::type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
