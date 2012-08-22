//=================================================================================================
/*!
//  \file blaze/math/shims/Reset.h
//  \brief Header file for the reset shim
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

#ifndef _BLAZE_MATH_SHIMS_RESET_H_
#define _BLAZE_MATH_SHIMS_RESET_H_


namespace blaze {

//=================================================================================================
//
//  RESET SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Resetting the given value/object to the default value.
// \ingroup math_shims
//
// \param resettable The value/object to be resetted.
// \return void
//
// The reset shim represents an abstract interface for the resetting of a value/object of
// any given data type to its default value. Values of built-in data type are reset to zero.
*/
template< typename Type >
inline void reset( Type& resettable )
{
   resettable = Type(0);
}
//*************************************************************************************************

} // namespace blaze

#endif
