//=================================================================================================
/*!
//  \file src/util/logging/Logger.cpp
//  \brief Source file for the Logger class
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


//*************************************************************************************************
// Platform/compiler-specific includes
//*************************************************************************************************

#include <blaze/system/WarningDisable.h>


//*************************************************************************************************
// MPI includes
//*************************************************************************************************

#include <blaze/system/MPI.h>
#if BLAZE_MPI_PARALLEL_MODE
#  include <mpi.h>
#endif


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstring>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <blaze/util/logging/Logger.h>


namespace blaze {

namespace logging {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Logger class.
*/
Logger::Logger()
   : Singleton<Logger,SystemClock>()  // Initialization of the Singleton base object
   , mutex_()                         // Synchronization mutex for thread-parallel logging
   , log_  ()                         // The log file
{}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the Logger class.
//
// The destructor of the Logger class writes the bottom line of the log file and closes
// the file.
*/
Logger::~Logger()
{
   if( log_.is_open() )
   {
      const std::time_t t = theSystemClock()->now();
      std::tm* localTime;
      char c[100];

      // Calculation of the local time
      localTime = std::localtime( &t );

      // Construction of the filename and time string
      std::strftime( c, 100, "%A, %d.%B %Y, %H:%M:%S", localTime );

      // Writing the bottom line
      const size_t length( std::strlen( c ) );
      log_ << "\n--LOG END, " << std::setw(length) << c
         << std::setw(90-length) << std::setfill('-') << "\n\n" << std::endl;

      // Closing the log file
      log_.close();
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Opens and initializes the log file.
//
// \return void
//
// This function is responsible for opening the log file and writing the header line. In case
// of a non-MPI-parallel simulation the function creates the log file 'blaze.log', which contains
// all logging information from all logging levels. In case of a MPI parallel simulation, each
// process creates his own individual log file called 'blazeX.log', where 'X' is replaced by
// the according rank the process has in the MPI_COMM_WORLD communicator.
*/
void Logger::openLogFile()
{
   // Creating the filename
   std::ostringstream filename;
   filename << "blaze";

#if BLAZE_MPI_PARALLEL_MODE
   int initialized( 0 );
   MPI_Initialized( &initialized );

   if( initialized ) {
      int rank( 0 );
      MPI_Comm_rank( MPI_COMM_WORLD, &rank );  // Estimating the rank of this process
      filename << rank;
   }
#endif

   filename << ".log";

   // Creating the time and date string for the header line
   const std::time_t t = theSystemClock()->start();
   std::tm* localTime;
   char c[100];

   localTime = std::localtime( &t );
   strftime( c, 100, "%A, %d.%B %Y, %H:%M:%S", localTime );

   // Opening the log file
   log_.open( filename.str().c_str(), std::ofstream::out | std::ofstream::trunc );
   if( !log_.is_open() ) {
      std::ostringstream oss;
      oss << " Error opening log file '" << filename << "' !\n";
      throw std::runtime_error( oss.str() );
   }

   // Writing the header line
   const size_t length( std::strlen( c ) );
   log_ << "\n\n--LOG BEGIN, " << std::setw(length) << c
        << std::setw(88-length) << std::setfill('-') << "\n\n" << std::endl;
}
//*************************************************************************************************

} // namespace logging

} // namespace blaze
