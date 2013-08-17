//=================================================================================================
/*!
//  \file blaze/util/mpl/Not.h
//  \brief Header file for the Not class template
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

#ifndef _BLAZE_UTIL_MPL_NOT_H_
#define _BLAZE_UTIL_MPL_NOT_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type negation.
// \ingroup mpl
//
// The Not class template negates the given compile time condition. In case the given condition
// would evaluate to \a true, the nested member enumeration is set to \a false and vice versa:

   \code
   using namespace blaze;

   Not< IsIntegral<int> >::value  // Evaluates to 0
   Not< IsDouble<int>   >::value  // Evaluates to 1
   \endcode
*/
template< typename C >  // Condition to be negated
struct Not
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   enum { value = !C::value };
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
