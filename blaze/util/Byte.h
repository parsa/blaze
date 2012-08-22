//=================================================================================================
/*!
//  \file blaze/util/Byte.h
//  \brief Header file for the byte type
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

#ifndef _BLAZE_UTIL_BYTE_H_
#define _BLAZE_UTIL_BYTE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Integral.h>
#include <blaze/util/constraints/Size.h>


namespace blaze {

//=================================================================================================
//
//  TYPE DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Byte data type of the Blaze library.
// \ingroup util
//
// The \a byte data type is guaranteed to be an integral data type of size 1.
*/
typedef unsigned char  byte;
//*************************************************************************************************




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace check {

// Constraint on the byte data type
BLAZE_CONSTRAINT_MUST_BE_INTEGRAL_TYPE( byte );
BLAZE_CONSTRAINT_MUST_HAVE_SIZE( byte, 1 );

} // namespace check
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
