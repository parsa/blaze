//=================================================================================================
/*!
//  \file blaze/util/SystemClock.h
//  \brief Header file for the SystemClock class
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

#ifndef _BLAZE_UTIL_SYSTEMCLOCK_H_
#define _BLAZE_UTIL_SYSTEMCLOCK_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <ctime>
#include <blaze/util/Null.h>
#include <blaze/util/singleton/Singleton.h>
#include <blaze/util/SystemClockID.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief System clock of the Blaze library.
// \ingroup util
//
// The SystemClock class represents the system clock of the Blaze library. The system clock
// is the central timing functionality that can be used to query for the start time of the
// process, the current timestamp and the elapsed time since the start of the process. The
// following example demonstrates how the single system clock instance is acquired via the
// theSystemClock() functcion and how the system clock can be used:

   \code
   // The single system clock instance is accessible via the theSystemClock() function
   SystemClockID systemClock = theSystemClock();

   time_t start   = systemClock->start();    // Querying the start time of the process
   time_t current = systemClock->current();  // Querying the current timestamp
   time_t elapsed = systemClock->elapsed();  // Querying the elapsed time
   \endcode
*/
class SystemClock : private Singleton<SystemClock>
{
 private:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SystemClock();
   //@}
   //**********************************************************************************************

 public:
   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~SystemClock();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline time_t start  () const;
   inline time_t now    () const;
   inline time_t elapsed() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static time_t start_;  //!< Timestamp for the start of the process.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   friend SystemClockID theSystemClock();
   BLAZE_BEFRIEND_SINGLETON;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  SYSTEM CLOCK SETUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name System clock setup functions */
//@{
inline SystemClockID theSystemClock();
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a handle to the Blaze system clock.
// \ingroup util
//
// \return Handle to the active system clock.
*/
inline SystemClockID theSystemClock()
{
   return SystemClock::instance();
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the timestamp for the start of the process.
//
// \return Timestamp for the start of the process.
*/
inline time_t SystemClock::start() const
{
   return start_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current timestamp.
//
// \return The current timestamp.
*/
inline time_t SystemClock::now() const
{
   return time( NULL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the elapsed time since the start of the process (in seconds).
//
// \return Elapsed time since the start of the process (in seconds).
*/
inline time_t SystemClock::elapsed() const
{
   return std::time( NULL ) - start_;
}
//*************************************************************************************************

} // namespace blaze

#endif
