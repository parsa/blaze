//=================================================================================================
/*!
//  \file blaze/util/logging/ProgressSection.h
//  \brief Header file for the log progress section
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

#ifndef _BLAZE_UTIL_LOGGING_PROGRESSSECTION_H_
#define _BLAZE_UTIL_LOGGING_PROGRESSSECTION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/logging/LogSection.h>


namespace blaze {

namespace logging {

//=================================================================================================
//
//  BLAZE_LOG_PROGRESS_SECTION MACRO
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Logging section for progress information.
// \ingroup logging
//
// This macro starts a log section for information messages. These messages are written to the
// log file(s) in case the blaze::loglevel has been set to \a progress or higher. The following
// example demonstrates how this log section is used:

   \code
   int main( int argc, char** argv )
   {
      // Initialization of the MPI system (for MPI parallel simulations)
      // The MPI system must be initialized before any logging functionality may be used. In
      // case it was not called before the first log section it is assumed that the simulation
      // does not run in parallel. Thus in MPI-parallel simulations it is strongly recommended
      // to make MPI_Init() the very first call of the main function.
      MPI_Init( &argc, &argv );

      // ...

      // Log section for progress information
      // This section is only executed in case the logging level is at least 'progress'. The
      // macro parameter specifies the name of the log handle (in this example 'log') that
      // can be used as a stream to log any kind of streamable information.
      BLAZE_LOG_PROGRESS_SECTION( log ) {
         log << " Only printed within an active BLAZE_LOG_PROGRESS_SECTION!\n";
      }

      // ...

      // Finalizing the MPI system (for MPI parallel simulations)
      // The MPI system must be finalized after the last pe functionality has been used. It
      // is recommended to make MPI_Finalize() the very last call of the main function.
      MPI_Finalize();
   }
   \endcode

// Note that uncaught exceptions emitted from the blaze::BLAZE_LOG_PROGRESS_SECTION might result in
// lost and/or unlogged information!
*/
#define BLAZE_LOG_PROGRESS_SECTION( NAME ) \
   if( blaze::logging::loglevel >= blaze::logging::progress ) \
      if( blaze::logging::LogSection NAME = blaze::logging::progress )
//*************************************************************************************************

} // namespace logging

} // namespace blaze

#endif
