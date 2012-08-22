//=================================================================================================
/*!
//  \file blaze/util/timing/Timing.h
//  \brief Timing module documentation
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

#ifndef _BLAZE_UTIL_TIMING_TIMING_H_
#define _BLAZE_UTIL_TIMING_TIMING_H_


namespace blaze {

//=================================================================================================
//
//  DOXYGEN DOCUMENTATION
//
//=================================================================================================

//*************************************************************************************************
//! Namespace for the time measurement module.
namespace timing {}
//*************************************************************************************************


//*************************************************************************************************
/*!\defgroup timing Time measurement
// \ingroup util
//
// \image html clock.png
// \image latex clock.eps "Timing submodule" width=200pt
//
// The timing submodule offers the necessary functionality for timing and benchmarking purposes.
// The central element of the timing module is the Timer class. Depending on a chosen timing
// policy, this class offers the possibility to measure both single times and time series. In
// order to make time measurement as easy as possible, the Blaze library offers the two classes
// WcTimer and CpuTimer (both using the Timer class) to measure both wall clock and CPU time.
// The following example gives an impression on how time measurement for a single time works
// with the the Blaze library. Note that in this example the WcTimer could be easily replaced
// with the CpuTimer if instead of the wall clock time the CPU time was to be measured.

   \code
   // Creating a new wall clock timer immediately starts a new time measurement
   WcTimer timer;

   ...  // Programm or code fragment to be measured

   // Stopping the time measurement
   timer.end();

   // Evaluation of the measured time
   double time = timer.last();
   \endcode

// As already mentioned, it is also possible to start several time measurements with a single
// timer to evaluate for instance the minimal, the maximal or the average time of a specific
// task. The next example demonstrates a possible setup for such a series of time measurements:

   \code
   // Creating a new wall clock timer
   WcTimer timer;

   ...  // Additional setup code

   // Starting 10 wall clock time measurements
   for( unsigned int i=0; i<10; ++i ) {
      timer.start();
      ...  // Programm or code fragment to be measured
      timer.end();
   }

   // After the measurements, the desired timing results can be calculated, as for instance the
   // average wall clock time
   double average = timer.average();
   \endcode
*/
//*************************************************************************************************

} // namespace blaze

#endif
