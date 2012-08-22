//=================================================================================================
/*!
//  \file blaze/util/timing/WcTimer.h
//  \brief Progress timer for wall clock time measurements
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

#ifndef _BLAZE_UTIL_TIMING_WCTIMER_H_
#define _BLAZE_UTIL_TIMING_WCTIMER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/timing/Timer.h>
#include <blaze/util/timing/WcPolicy.h>


namespace blaze {

namespace timing {

//=================================================================================================
//
//  TYPE DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Progress timer for wall clock time measurements.
// \ingroup timing
//
// The WcTimer combines the Timer class template with the WcPolicy timing policy. It measures
// the amount of "wall clock" time elapsing for the processing of a programm or code fragment.
// In contrast to the measurement of CPU time, the wall clock time also contains waiting times
// such as input/output operations.
*/
typedef Timer<WcPolicy>  WcTimer;
//*************************************************************************************************

} // timing

} // blaze

#endif
