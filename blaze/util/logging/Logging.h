//=================================================================================================
/*!
//  \file blaze/util/logging/Logging.h
//  \brief Logging module documentation
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

#ifndef _BLAZE_UTIL_LOGGING_LOGGING_H_
#define _BLAZE_UTIL_LOGGING_LOGGING_H_


namespace blaze {

//=================================================================================================
//
//  DOXYGEN DOCUMENTATION
//
//=================================================================================================

//*************************************************************************************************
//! Namespace for the logging module.
namespace logging {}
//*************************************************************************************************


//*************************************************************************************************
/*!\defgroup logging Logging
// \ingroup util
//
// The logging submodule offers the functionality for the creation of log information in both
// non-parallel and MPI-/thread-parallel environments. The logging functionality is implemented
// such that in case no logging is required no runtime and memory overhead occur. However, in
// case it is necessary to log information, this is done as efficiently and reliably as possible.
// In non-parallel environments, a single log file named 'blaze.log' is created, containing all
// the information of the single process. In MPI-parallel environments, each process creates
// his own log file named 'blazeX.log', where 'X' is replaced by his according process rank in
// the MPI_COMM_WORLD communicator. Depending on the selected logging level information about
// (severe) errors, warnings, important information, progress reports, debug information, and
// detailed output is written to the log file(s). The global log level is specified via the
// blaze::loglevel variable in the configuration file "blaze/config/Logging.h". The following
// logging level's are available:
//
//  - \a inactive: Completely deactivates the logging functionality, i.e., no log file(s) will
//                 be written. Since this setting can immensely complicate error correction, it
//                 is not recommended to use this setting!
//  - \a error   : Only (severe) errors are written to the log file(s).
//  - \a warning : Extends the \a error setting by warning messages.
//  - \a info    : Extends the \a warning setting by additional informative messages (default).
//  - \a progress: Extends the \a info setting by progress information.
//  - \a debug   : Extends the \a progress setting by debug information.
//  - \a detail  : Extends the \a debug setting by very fine grained detail information.
//
// Logging is done via one of the six following log sections:
//
//  - blaze::BLAZE_LOG_ERROR_SECTION   : Section for (severe) error messages; blaze::loglevel >= \a error
//  - blaze::BLAZE_LOG_WARNING_SECTION : Section for warning messages; blaze::loglevel >= \a warning
//  - blaze::BLAZE_LOG_INFO_SECTION    : Section for important information; blaze::loglevel >= \a info
//  - blaze::BLAZE_LOG_PROGRESS_SECTION: Section for progress information; blaze::loglevel >= \a progress
//  - blaze::BLAZE_LOG_DEBUG_SECTION   : Section for debug information; blaze::loglevel >= \a debug
//  - blaze::BLAZE_LOG_DETAIL_SECTION  : Section for detail information; blaze::loglevel >= \a detail
//
// The following example demonstrates the use of the log sections:

   \code
   int main( int argc, char** argv )
   {
      // Initialization of the MPI system (for MPI parallel simulations only)
      // The MPI system must be initialized before any logging functionality may be used. In
      // case it was not called before the first log section it is assumed that the simulation
      // does not run in parallel. Thus in MPI-parallel simulations it is strongly recommended
      // to make MPI_Init() the very first call of the main function.
      MPI_Init( &argc, &argv );

      // ...

      // Log section for error messages
      // This section is only executed in case the logging level is at least 'error'. The
      // macro parameter specifies the name of the log handle (in this example 'log') that
      // can be used as a stream to log any kind of streamable information.
      BLAZE_LOG_ERROR_SECTION( log ) {
         log << " Only printed within an active BLAZE_LOG_ERROR_SECTION\n"
             << "   for demonstration purposes!\n";
      }

      // Log section for warning messages
      // This section is only executed in case the logging level is at least 'warning'. Again,
      // the macro parameter specifies the name of the log handle that can be used as a stream
      // to any kind of streamable information.
      BLAZE_LOG_WARNING_SECTION( log ) {
         log << " Only printed within an active BLAZE_LOG_WARNING_SECTION!\n";
      }

      // ...

      // Finalizing the MPI system (for MPI parallel simulations only)
      // The MPI system must be finalized after the last pe functionality has been used. It
      // is recommended to make MPI_Finalize() the very last call of the main function.
      MPI_Finalize();
   }
   \endcode

// Executing this main() function results in the following log file (provided the logging level
// is set accordingly):

   \code
   --LOG BEGIN, Thursday, 14.January 2010, 08:31--------------------------------------------


   [ERROR   ]  Only printed within an active BLAZE_LOG_ERROR_SECTION
                 for demonstration purposes!

   [WARNING ]  Only printed within an active BLAZE_LOG_WARNING_SECTION!


   --LOG END, Thursday, 14.January 2010, 08:31----------------------------------------------
   \endcode

// The next example shows to nested log sections. Please note that this combination only makes
// sense in case the outer log section has a lower log level:

   \code
   // ...

   BLAZE_LOG_ERROR_SECTION( error )
   {
      error << " Only printed within an active BLAZE_LOG_ERROR_SECTION\n";

      BLAZE_LOG_WARNING_SECTION( warning ) {
         warning << " Only printed within an active BLAZE_LOG_WARNING_SECTION\n";
      }

      error << " Again only printed within an active BLAZE_LOG_ERROR_SECTION\n";
   }

   // ...
   \endcode

// Although the BLAZE_LOG_ERROR_SECTION is the outer section, the BLAZE_LOG_WARNING_SECTION
// appears first in the log file:

   \code
   --LOG BEGIN, Thursday, 14.January 2010, 08:31--------------------------------------------


   [WARNING ]  Only printed within an active BLAZE_LOG_WARNING_SECTION

   [ERROR   ]  Only printed within an active BLAZE_LOG_ERROR_SECTION
               Again only printed within an active BLAZE_LOG_ERROR_SECTION


   --LOG END, Thursday, 14.January 2010, 08:31----------------------------------------------
   \endcode

// In order to commit the messages in the correct order, it is necessary to manually call the
// commit function:

   \code
   // ...

   BLAZE_LOG_ERROR_SECTION( error )
   {
      error << " Only printed within an active BLAZE_LOG_ERROR_SECTION\n";

      // Manual call of the commit function
      error.commit();

      BLAZE_LOG_WARNING_SECTION( warning ) {
         warning << " Only printed within an active BLAZE_LOG_WARNING_SECTION\n";
      }

      error << " Again only printed within an active BLAZE_LOG_ERROR_SECTION\n";
   }

   // ...
   \endcode

// This results in the following log file:

   \code
   --LOG BEGIN, Thursday, 14.January 2010, 08:31--------------------------------------------


   [ERROR   ]  Only printed within an active BLAZE_LOG_ERROR_SECTION

   [WARNING ]  Only printed within an active BLAZE_LOG_WARNING_SECTION

   [ERROR   ]  Again only printed within an active BLAZE_LOG_ERROR_SECTION


   --LOG END, Thursday, 14.January 2010, 08:31----------------------------------------------
   \endcode
*/
//*************************************************************************************************

} // namespace blaze

#endif
