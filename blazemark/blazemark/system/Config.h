//=================================================================================================
/*!
//  \file blazemark/system/Config.h
//  \brief General settings for the blaze benchmark suite
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

#ifndef _BLAZEMARK_SYSTEM_CONFIG_H_
#define _BLAZEMARK_SYSTEM_CONFIG_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/StaticAssert.h>


namespace blazemark {

//=================================================================================================
//
//  GENERAL CONFIGURATION
//
//=================================================================================================

#include <blazemark/config/Config.h>




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( reps > 0 );
BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( element_t );

}
/*! \endcond */
//*************************************************************************************************

} // namespace blazemark

#endif
