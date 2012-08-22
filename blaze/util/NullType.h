//=================================================================================================
/*!
//  \file blaze/util/NullType.h
//  \brief Utility type for generic codes
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

#ifndef _BLAZE_UTIL_NULLTYPE_H_
#define _BLAZE_UTIL_NULLTYPE_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Utility type for generic codes.
// \ingroup util
//
// The NullType class represents an invalid or terminating data type for generic codes. For
// instance, the TypeList class uses the NullType as terminating data type for the type list.
*/
class NullType
{};
//*************************************************************************************************

} // namespace blaze

#endif
