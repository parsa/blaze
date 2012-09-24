//=================================================================================================
/*!
//  \file blaze/math/shims/Square.h
//  \brief Header file for the square shim
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

#ifndef _BLAZE_MATH_SHIMS_SQUARE_H_
#define _BLAZE_MATH_SHIMS_SQUARE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/traits/MultExprTrait.h>


namespace blaze {

//=================================================================================================
//
//  SQUARE SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Squaring the given value/object.
// \ingroup math_shims
//
// \param a The value/object to be squared.
// \return The result of the square operation.
//
// The square shim represents an abstract interface for squaring a value/object of any given
// data type. For values of built-in data type this results in a plain multiplication.
*/
template< typename T >
inline const typename MultExprTrait<T,T>::Type sq( const T& a )
{
   return a * a;
}
//*************************************************************************************************

} // namespace blaze

#endif
