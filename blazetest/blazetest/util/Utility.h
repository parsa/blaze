//=================================================================================================
/*!
//  \file blazetest/util/Utility.h
//  \brief Header file for utility test functionality
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

#ifndef _BLAZETEST_UTIL_UTILITY_H_
#define _BLAZETEST_UTIL_UTILITY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/shims/Equal.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blazetest {

//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary comparison function.
//
// \param m1 The left-hand side element to be compared.
// \param m2 The right-hand side element to be compared.
// \return \a true in case the two elements are equal, \a false in case the elements are not equal.
//
// This function performs an excessive comparison of the two given elements. It utilizes both
// equality and inequality comparison and additionally swaps the two elements.
*/
template< typename T1    // Type of the first element
        , typename T2 >  // Type of the second element
bool isEqual( const T1& m1, const T2& m2 )
{
   if( blaze::IsNumeric<T1>::value && blaze::IsNumeric<T2>::value )
      return blaze::equal( m1, m2 );
   else
      return ( m1 == m2 && !( m1 != m2 ) && m2 == m1 && !( m2 != m1 ) );
}
//*************************************************************************************************

} // namespace blazetest

#endif
