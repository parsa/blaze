//=================================================================================================
/*!
//  \file blaze/util/timing/WcPolicy.h
//  \brief Wall clock timing policy for the Timer class.
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

#ifndef _BLAZE_UTIL_TIMING_WCPOLICY_H_
#define _BLAZE_UTIL_TIMING_WCPOLICY_H_


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
/*!\brief Timing policy for the measurement of the wall clock time.
// \ingroup timing
//
// The WcPolicy class represents the timing policy for wall clock time measurements that can be
// used in combination with the Timer class template. This combination is realized with the WcTimer
// type definition.
*/
struct WcPolicy
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
/*!\brief Returns a timestamp of the current wall clock time in seconds.
//
// \return Wall clock timestamp in seconds.
*/
inline double WcPolicy::getTimestamp()
{
   return getWcTime();
}
//*************************************************************************************************

} // timing

} // blaze

#endif
