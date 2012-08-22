//=================================================================================================
/*!
//  \file blaze/math/shims/Clear.h
//  \brief Header file for the clear shim
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

#ifndef _BLAZE_MATH_SHIMS_CLEAR_H_
#define _BLAZE_MATH_SHIMS_CLEAR_H_


namespace blaze {

//=================================================================================================
//
//  CLEAR SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Clearing the given value/object to the default state.
// \ingroup math_shims
//
// \param clearable The value/object to be cleared.
// \return void
//
// The clear shim represents an abstract interface for clearing a value/object of any given
// data type to its default state. Values of built-in data type are reset to zero.
*/
template< typename Type >
inline void clear( Type& clearable )
{
   clearable = Type(0);
}
//*************************************************************************************************

} // namespace blaze

#endif
