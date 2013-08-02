//=================================================================================================
/*!
//  \file blaze/util/AlignmentCheck.h
//  \brief Header file for the alignment check function
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

#ifndef _BLAZE_UTIL_ALIGNMENTCHECK_H_
#define _BLAZE_UTIL_ALIGNMENTCHECK_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/AlignmentTrait.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  SIZETRAIT CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks the alignment of the given
// \ingroup util
//
// \param address The address to be checked.
// \return \a true in case the address is properly aligned, \a false if it is not.
//
// This function performs an alignment check on the given address. For instance, for fundamental
// data types that can be vectorized via SSE or AVX instructions, the proper alignment is 16 or
// 32 bytes, respectively. In case the given address is properly aligned, the function returns
// \a true, otherwise it returns \a false.
*/
template< typename T >
bool checkAlignment( const T* address )
{
   return !( reinterpret_cast<size_t>( address ) % AlignmentTrait<T>::value );
}
//*************************************************************************************************

} // namespace blaze

#endif
