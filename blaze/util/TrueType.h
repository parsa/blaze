//=================================================================================================
/*!
//  \file blaze/util/TrueType.h
//  \brief Header file for the TrueType type/value trait base class
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

#ifndef _BLAZE_UTIL_TRUETYPE_H_
#define _BLAZE_UTIL_TRUETYPE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/type_traits/integral_constant.hpp>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type traits base class.
// \ingroup util
//
// The TrueType class is used as base class for type traits and value traits that evaluate to
// \a true.
*/
typedef boost::true_type  TrueType;
//*************************************************************************************************

} // namespace blaze

#endif
