//=================================================================================================
/*!
//  \file blaze/util/typetraits/AddVolatile.h
//  \brief Header file for the AddVolatile type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ADDVOLATILE_H_
#define _BLAZE_UTIL_TYPETRAITS_ADDVOLATILE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/add_volatile.hpp>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Addition of a top level 'volatile' qualifier.
// \ingroup type_traits
//
// The AddVolatile type trait adds a top level 'volatile' qualifier to the given type \a T.

   \code
   blaze::AddVolatile<int>::Type           // Results in 'int volatile'
   blaze::AddVolatile<int*>::Type          // Results in 'int* volatile'
   blaze::AddVolatile<int&>::Type          // Results in 'int&'
   blaze::AddVolatile<int volatile>::Type  // Results in 'int volatile'
   blaze::AddVolatile<int const>::Type     // Results in 'int const volatile'
   \endcode
*/
template< typename T >
struct AddVolatile
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   typedef typename boost::add_volatile<T>::type  Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
