//=================================================================================================
/*!
//  \file src/util/logging/Logger.cpp
//  \brief Source file for the Logger class
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
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
