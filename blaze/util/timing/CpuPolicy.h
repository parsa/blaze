//=================================================================================================
/*!
//  \file blaze/util/timing/CpuPolicy.h
//  \brief CPU timing policy for the Timer class.
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

#ifndef _BLAZE_UTIL_TIMING_CPUPOLICY_H_
#define _BLAZE_UTIL_TIMING_CPUPOLICY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Time.h>


namespace blaze {

namespace timing {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Timing policy for the measurement of the CPU time.
// \ingroup timing
//
// The CpuPolicy class represents the timing policy for CPU time measurements that can be used
// in combination with the Timer class template. This combination is realized with the CpuTimer
// type definition.
*/
struct CpuPolicy
{
 public:
   //**Timing functions****************************************************************************
   /*!\name Timing functions */
   //@{
   static inline double getTimestamp();
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TIMING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a timestamp of the current CPU time in seconds.
//
// \return CPU timestamp in seconds.
*/
inline double CpuPolicy::getTimestamp()
{
   return getCpuTime();
}
//*************************************************************************************************

} // timing

} // blaze

#endif
