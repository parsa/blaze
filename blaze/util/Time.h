//=================================================================================================
/*!
//  \file blaze/util/Time.h
//  \brief Header file for time functions
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

#ifndef _BLAZE_UTIL_TIME_H_
#define _BLAZE_UTIL_TIME_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#if defined(_MSC_VER)
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#  include <winsock.h>
#  include <time.h>
#  include <sys/timeb.h>
#else
#  include <sys/resource.h>
#  include <sys/time.h>
#  include <sys/types.h>
#endif
#include <ctime>
#include <string>


namespace blaze {

//=================================================================================================
//
//  TIME FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Time functions */
//@{
inline std::string getDate();
inline std::string getTime();
inline double      getWcTime();
inline double      getCpuTime();
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a formated date string in the form YYYY-MM-DD
// \ingroup util
//
// \return Formated date string
*/
inline std::string getDate()
{
   std::time_t t;
   std::tm* localTime;
   char c[50];

   std::time( &t );
   localTime = std::localtime( &t );
   std::strftime( c, 50, "%Y-%m-%d", localTime );

   return std::string( c );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a formated time and date string
// \ingroup util
//
// \return Formated time and date string in the format WEEKDAY DAY.MONTH YEAR, HOUR:MINUTES
*/
inline std::string getTime()
{
   std::time_t t;
   std::tm* localTime;
   char c[50];

   std::time( &t );
   localTime = std::localtime( &t );
   std::strftime( c, 50, "%A, %d.%B %Y, %H:%M", localTime );

   return std::string( c );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current wall clock time in seconds.
// \ingroup util
//
// \return The current wall clock time in seconds.
*/
inline double getWcTime()
{
#ifdef WIN32
   struct _timeb timeptr;
   _ftime( &timeptr );
   return ( static_cast<double>( timeptr.time ) + static_cast<double>( timeptr.millitm )/1E3 );
#else
   struct timeval tp;
   gettimeofday( &tp, NULL );
   return ( static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6 );
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current CPU time in seconds.
// \ingroup util
//
// \return The current CPU time in seconds.
*/
inline double getCpuTime()
{
#ifdef WIN32
   FILETIME CreateTime, ExitTime, KernelTime, UserTime;
   SYSTEMTIME SysTime;

   if( GetProcessTimes( GetCurrentProcess(), &CreateTime, &ExitTime, &KernelTime, &UserTime ) != TRUE ) {
      return 0.0;
   }
   else {
      FileTimeToSystemTime( &UserTime, &SysTime );
      return ( static_cast<double>( SysTime.wSecond ) + static_cast<double>( SysTime.wMilliseconds )/1E3 );
   }
#else
   struct rusage ruse;
   getrusage( RUSAGE_SELF, &ruse );
   return ( static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6 );
#endif
}
//*************************************************************************************************

} // namespace blaze

#endif
