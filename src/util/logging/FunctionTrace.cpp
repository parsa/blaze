//=================================================================================================
/*!
//  \file src/util/logging/FunctionTrace.cpp
//  \brief Source file for the FunctionTrace class
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
// Includes
//*************************************************************************************************

#include <sstream>
#include <boost/format.hpp>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/logging/Logger.h>
#include <blaze/util/SystemClock.h>


namespace blaze {

namespace logging {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the FunctionTrace class.
//
// \param file The name of the file the traced function is contained in
// \param function The name of the traced function
*/
FunctionTrace::FunctionTrace( const std::string& file, const std::string& function )
   : file_    ( file     )  // The file name the traced function is contained in
   , function_( function )  // The name of the traced function
{
   // Writing the trace information
   std::ostringstream oss;
   oss << "[TRACE   ]";

   // Writing the elapsed time
   const time_t elapsed( theSystemClock()->elapsed() );
   const time_t hours  ( elapsed / time_t( 3600 ) );
   const time_t minutes( ( elapsed % time_t( 3600 ) ) / time_t( 60 ) );
   const time_t seconds( elapsed % time_t( 60 ) );
   oss << boost::format( "[%03d:%02d:%02d]" ) % hours % minutes % seconds;

   // Writing the message string
   oss << " + Entering function '" << function_ << "' in file '" << file_ << "'\n";

   // Logging the message string
   boost::shared_ptr<Logger> logger( Logger::instance() );
   logger->log( oss.str() );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the FunctionTrace class.
*/
FunctionTrace::~FunctionTrace()
{
   // Writing the trace information
   std::ostringstream oss;
   oss << "[TRACE   ]";

   // Writing the elapsed time
   const time_t elapsed( theSystemClock()->elapsed() );
   const time_t hours  ( elapsed / time_t( 3600 ) );
   const time_t minutes( ( elapsed % time_t( 3600 ) ) / time_t( 60 ) );
   const time_t seconds( elapsed % time_t( 60 ) );
   oss << boost::format( "[%03d:%02d:%02d]" ) % hours % minutes % seconds;

   // Writing the message string
   oss << " - Leaving function '" << function_ << "' in file '" << file_ << "'\n";

   // Logging the message string
   boost::shared_ptr<Logger> logger( Logger::instance() );
   logger->log( oss.str() );
}
//*************************************************************************************************

} // namespace logging

} // namespace blaze
