//=================================================================================================
/*!
//  \file blaze/math/shims/IsDefault.h
//  \brief Header file for the isDefault shim
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

#ifndef _BLAZE_MATH_SHIMS_ISDEFAULT_H_
#define _BLAZE_MATH_SHIMS_ISDEFAULT_H_


namespace blaze {

//=================================================================================================
//
//  ISDEFAULT SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the given value/object is in default state.
// \ingroup math_shims
//
// \param v The value/object to be tested for its default state.
// \return \a true in case the given value/object is in its default state, \a false otherwise.
//
// The isDefault shim represents an abstract interface for testing a value/object whether
// it is in its default state or not. In case the value/object is in its default state, the
// function returns \a true, otherwise it returns \a false. For built-in data types, the
// function returns \a true in case the current value is zero.

   \code
   const int i = 0;          // isDefault( i ) returns true
   double    d = 2.0;        // isDefault( d ) returns false
   Vec3      v1;             // isDefault( v1 ) returns true
   Vec3      v2( 0, 0, 0 );  // isDefault( v2 ) returns true since (0,0,0) is the default state
   Vec3      v3( 1, 2, 3 );  // isDefault( v3 ) returns false
   \endcode
*/
template< typename Type >
inline bool isDefault( const Type& v )
{
   return v == Type(0);
}
//*************************************************************************************************

} // namespace blaze

#endif
