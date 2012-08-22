//=================================================================================================
/*!
//  \file blaze/system/Solvers.h
//  \brief System settings for for the mathematical solvers
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

#ifndef _BLAZE_SYSTEM_SOLVER_H_
#define _BLAZE_SYSTEM_SOLVER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Precision.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>




//=================================================================================================
//
//  MATHEMATICAL SOLVER SETTINGS
//
//=================================================================================================

#include <blaze/config/Solvers.h>




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( blaze::solvers::maxIterations > 0 );

}
/*! \endcond */
//*************************************************************************************************

#endif
