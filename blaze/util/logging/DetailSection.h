//=================================================================================================
/*!
//  \file blaze/util/logging/DetailSection.h
//  \brief Header file for the log detail section
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

#ifndef _BLAZE_UTIL_LOGGING_DETAILSECTION_H_
#define _BLAZE_UTIL_LOGGING_DETAILSECTION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/logging/LogSection.h>


namespace blaze {

namespace logging {

//=================================================================================================
//
//  BLAZE_LOG_DETAIL_SECTION MACRO
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Logging section for debug information.
// \ingroup logging
//
// This macro starts a log section for detail information. These messages are written to the
// log file(s) in case the blaze::loglevel has been set to \a detail or higher. The following
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

      // Log section for detail information
      // This section is only executed in case the logging level is at least 'detail'. The
      // macro parameter specifies the name of the log handle (in this example 'log') that
      // can be used as a stream to log any kind of streamable information.
      BLAZE_LOG_DETAIL_SECTION( log ) {
         log << " Only printed within an active BLAZE_LOG_DETAIL_SECTION!\n";
      }

      // ...

      // Finalizing the MPI system (for MPI parallel simulations)
      // The MPI system must be finalized after the last pe functionality has been used. It
      // is recommended to make MPI_Finalize() the very last call of the main function.
      MPI_Finalize();
   }
   \endcode

// Note that uncaught exceptions emitted from the blaze::BLAZE_LOG_DETAIL_SECTION might result
// in lost and/or unlogged information!
*/
#define BLAZE_LOG_DETAIL_SECTION( NAME ) \
   if( blaze::logging::loglevel >= blaze::logging::detail ) \
      if( blaze::logging::LogSection NAME = blaze::logging::detail )
//*************************************************************************************************

} // namespace logging

} // namespace blaze

#endif
