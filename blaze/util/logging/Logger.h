//=================================================================================================
/*!
//  \file blaze/util/logging/Logger.h
//  \brief Header file for the Logger class
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

#ifndef _BLAZE_UTIL_LOGGING_LOGGER_H_
#define _BLAZE_UTIL_LOGGING_LOGGER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <fstream>
#include <boost/thread/mutex.hpp>
#include <blaze/util/singleton/Singleton.h>
#include <blaze/util/SystemClock.h>


namespace blaze {

namespace logging {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of a logger class.
// \ingroup logging
//
// The Logger class represents the core of the logging functionality. It is responsible for
// commiting logging messages immediately to the according log file(s). The logger works for
// both serial as well as MPI parallel environments. In case of a non-MPI-parallel simulation
// the Logger creates the log file 'blaze.log', which contains all logging information from all
// logging levels. In case of a MPI parallel simulation, each process creates his own individual
// log file called 'blazeX.log', where 'X' is replaced by the according rank the process has in
// the MPI_COMM_WORLD communicator.\n
// Note that the log file(s) are only created in case any logging information is created. This
// might for instance result in only a small number of log file(s) in MPI parallel simulations
// when only some of the processes encounter errors/warnings/etc.\n
// Note that the logging functionality may not be used before MPI_Init() has been finished. In
// consequence, this means that no global data that is initialized before the main() function
// may contain any use of the logging functionality!
*/
class Logger : private Singleton<Logger,SystemClock>
{
 private:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Logger();
   //@}
   //**********************************************************************************************

 public:
   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Logger();
   //@}
   //**********************************************************************************************

 private:
   //**Logging functions***************************************************************************
   /*!\name Logging functions */
   //@{
   template< typename Type > void log( const Type& message );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void openLogFile();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   boost::mutex  mutex_;  //!< Synchronization mutex for thread-parallel logging.
   std::ofstream log_;    //!< The log file.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   friend class LogSection;
   BLAZE_BEFRIEND_SINGLETON;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  LOGGING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Writes the log message to the log file.
//
// \param message The log message to be logged.
// \return void
//
// This function immediately commits the log message to the log file. The first call to this
// function will create the log file.
*/
template< typename Type >  // Type of the log message
void Logger::log( const Type& message )
{
   boost::mutex::scoped_lock lock( mutex_ );
   if( !log_.is_open() )
      openLogFile();
   log_ << message;
   log_.flush();
}
//*************************************************************************************************

} // namespace logging

} // namespace blaze

#endif
