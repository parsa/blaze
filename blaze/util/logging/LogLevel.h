//=================================================================================================
/*!
//  \file blaze/util/logging/LogLevel.h
//  \brief Header file for the logging levels
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

#ifndef _BLAZE_UTIL_LOGGING_LOGLEVEL_H_
#define _BLAZE_UTIL_LOGGING_LOGLEVEL_H_


namespace blaze {

namespace logging {

//=================================================================================================
//
//  LOGGING LEVELS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Logging levels.
// \ingroup logging
//
// The LogLevel type enumeration represents the type of the global logging level. It defines
// all possible levels for the logging functionality. Depending on the setting of the global
// logging level (see blaze::logLevel), more or less information will be written to the log
// file(s). The following logging levels are available:
//
//  - \a inactive: Completely deactivates the logging functionality, i.e., no log file(s) will
//                 be written. Since this setting can immensely complicate error correction, it
//                 is not recommended to use this setting!
//  - \a error   : Only (severe) errors are written to the log file(s).
//  - \a warning : Extends the \a error setting by warning messages (default).
//  - \a info    : Extends the \a warning setting by additional informative messages.
//  - \a progress: Extends the \a info setting by progress information.
//  - \a debug   : Extends the \a progress setting by debug information.
//  - \a detail  : Extends the \a debug setting by very fine grained detail information.
//
// \a inactive plays a special role in the way that it switches off the logging functionality
// completely, i.e., no log file(s) will be created. The highest logging level is \a error,
// which exclusively writes severe errors to the log file(s). The lowest logging level is
// \a detail, which can create a tremendous amount of logging information. Note that each
// logging level comprises all higher logging levels. For instance, \a progress will also
// print all errors and warning to the log file(s).
*/
enum LogLevel
{
   inactive = 0,  //!< Log level for no logging.
   error    = 1,  //!< Log level for (sever) errors.
   warning  = 2,  //!< Log level for warnings.
   info     = 3,  //!< Log level for high-level information.
   progress = 4,  //!< Log level for progress information.
   debug    = 5,  //!< Log level for debug information.
   detail   = 6   //!< Log level for detail information.
};
//*************************************************************************************************

} // namespace logging

} // namespace blaze

#endif
